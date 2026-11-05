"""BindingDB API requests, validation, and retry policy."""

import json
import logging
import re
from datetime import timezone
from email.utils import parsedate_to_datetime
from importlib.metadata import version

import httpx

from ligandscout.models import Measurement, SearchOutcome, SourceError, SourceResult

URL = "https://www.bindingdb.org/rest/getLigandsByUniprots"
LOG = logging.getLogger(__name__)


class SchemaError(ValueError):
    """The source response does not match the BindingDB contract."""


class EmptyResponseError(SchemaError):
    """An empty response cannot establish whether records exist."""


def reject_constant(value):
    raise SchemaError(f"Invalid JSON constant: {value}")


def strict_json(body: bytes):
    return json.loads(body, parse_constant=reject_constant)


def record_key(record: Measurement):
    return tuple(
        value or ""
        for value in (
            record.compound_id,
            record.measurement_type,
            record.measurement_value,
            record.pmid,
            record.doi,
        )
    )


def parse_records(body: bytes) -> list[Measurement]:
    if not body.strip():
        raise EmptyResponseError("BindingDB returned an ambiguous empty response.")
    data = strict_json(body)
    if data == "":
        raise EmptyResponseError("BindingDB returned an ambiguous empty response.")
    if not isinstance(data, dict):
        raise SchemaError("Expected a JSON object.")
    envelope = data.get("getLindsByUniprotsResponse")
    if not isinstance(envelope, dict) or not isinstance(
        envelope.get("affinities"), list
    ):
        raise SchemaError("Expected getLindsByUniprotsResponse.affinities list.")
    records = []
    for index, row in enumerate(envelope["affinities"]):
        if not isinstance(row, dict):
            raise SchemaError(f"Record {index} must be an object.")
        for field in ("monomerid", "affinity_type", "affinity"):
            if not isinstance(row.get(field), str) or not row[field].strip():
                raise SchemaError(f"Record {index}: {field} must be a nonempty string.")
        if not re.fullmatch(r"[0-9]+", row["monomerid"]):
            raise SchemaError(f"Record {index}: monomerid must contain ASCII digits.")
        optional = {}
        for field in ("query", "smile", "pmid", "doi"):
            value = row.get(field)
            if value is not None and not isinstance(value, str):
                raise SchemaError(f"Record {index}: {field} must be a string or null.")
            optional[field] = value if value is not None and value.strip() else None
        records.append(
            Measurement(
                compound_id=row["monomerid"],
                compound_name=None,
                smiles=optional["smile"],
                inchikey=None,
                target_description=optional["query"],
                evidence_type="reported_activity",
                measurement_type=row["affinity_type"],
                measurement_value=row["affinity"],
                measurement_unit=None,
                pmid=optional["pmid"],
                doi=optional["doi"],
                source_url="https://www.bindingdb.org/rwd/bind/chemsearch/marvin/MolStructure.jsp?monomerid="
                + row["monomerid"],
            )
        )
    return sorted(records, key=record_key)


def retry_delay(response, attempt, now):
    delay = float(2**attempt)
    value = response.headers.get("Retry-After", "")
    if re.fullmatch(r"[0-9]+", value):
        # Avoid float overflow on unreasonable server values.
        digits = value.lstrip("0") or "0"
        requested = 61 if len(digits) > 2 else min(int(digits), 61)
    else:
        try:
            date = parsedate_to_datetime(value)
            if date.tzinfo is None:
                date = date.replace(tzinfo=timezone.utc)
            requested = max(0, (date - now()).total_seconds())
        except (ValueError, TypeError, OverflowError):
            requested = 0
    return max(delay, requested)


def request_once(client, accession):
    response = client.get(
        URL,
        params={"uniprot": accession, "response": "application/json"},
        headers={
            "Accept": "application/json",
            "User-Agent": f"ligandscout/{version('ligandscout')}",
        },
        timeout=30.0,
        follow_redirects=False,
    )
    for hop in range(6):
        if response.status_code not in (301, 302, 303, 307, 308):
            return response, None
        location = response.headers.get("Location")
        try:
            destination = response.url.join(location) if location else None
        except (httpx.InvalidURL, ValueError):
            destination = None
        if (
            hop == 5
            or destination is None
            or destination.scheme != "https"
            or not destination.host
            or destination.userinfo
        ):
            return (
                response,
                "BindingDB redirect policy failed (invalid destination or redirect limit).",
            )
        response = client.get(
            destination,
            headers={
                "Accept": "application/json",
                "User-Agent": f"ligandscout/{version('ligandscout')}",
            },
            timeout=30.0,
            follow_redirects=False,
        )
    raise AssertionError("Unreachable redirect state")


def fetch_bindingdb(accession, *, client, sleep, now) -> SearchOutcome:
    response = None
    error = None
    records = []
    LOG.info("Searching BindingDB for %s", accession)
    for attempt in range(3):
        response = None
        try:
            response, redirect_error = request_once(client, accession)
        except httpx.RequestError as exc:
            error = SourceError(
                "transport", f"BindingDB request failed ({type(exc).__name__}).", None
            )
            if (
                isinstance(exc, (httpx.TimeoutException, httpx.NetworkError))
                and attempt < 2
            ):
                delay = 2**attempt
            else:
                break
        else:
            if redirect_error:
                error = SourceError("http", redirect_error, response.status_code)
                break
            if response.status_code == 200:
                try:
                    records = parse_records(response.content)
                except EmptyResponseError:
                    error = SourceError(
                        "schema",
                        "BindingDB returned an empty response; cannot determine whether matches exist. Retry later.",
                        200,
                    )
                    if attempt == 2:
                        break
                    LOG.warning(
                        "BindingDB %s returned an ambiguous empty HTTP 200 response",
                        accession,
                    )
                    delay = retry_delay(response, attempt, now)
                    if delay > 60:
                        break
                except (ValueError, UnicodeError, RecursionError):
                    LOG.debug("BindingDB parsing failed", exc_info=True)
                    error = SourceError(
                        "schema",
                        "BindingDB response has invalid JSON or record schema.",
                        200,
                    )
                    break
                else:
                    error = None
                    break
            else:
                error = SourceError(
                    "http",
                    f"BindingDB returned HTTP {response.status_code}.",
                    response.status_code,
                )
                if (
                    response.status_code != 429
                    and not 500 <= response.status_code <= 599
                ) or attempt == 2:
                    break
                delay = retry_delay(response, attempt, now)
                if delay > 60:
                    error.message += " Retry later: requested wait exceeds 60 seconds."
                    break
        LOG.warning(
            "BindingDB %s: retry %s/3 in %s seconds", accession, attempt + 2, delay
        )
        sleep(delay)
    body = response.content if response is not None else None
    result = SourceResult(
        schema_version=1,
        accession=accession,
        source="BindingDB",
        retrieved_at=now()
        .astimezone(timezone.utc)
        .isoformat(timespec="microseconds")
        .replace("+00:00", "Z"),
        request={
            "url": URL,
            "parameters": {"uniprot": accession, "response": "application/json"},
        },
        source_version=None,
        status="failed" if error else ("success" if records else "no_matches"),
        records=records,
        error=error,
    )
    LOG.info(
        "BindingDB %s: %s; %s compounds, %s measurements",
        accession,
        result.status,
        len({r.compound_id for r in records}),
        len(records),
    )
    return SearchOutcome(result, body)
