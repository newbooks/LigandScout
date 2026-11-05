"""Audited, throttled JSON requests for annotation and chemical enrichment."""

import logging
from dataclasses import dataclass, field
from importlib.metadata import version

import httpx

from ligandscout.bindingdb import retry_delay, strict_json
from ligandscout.models import HttpExchange, SourceError

LOG = logging.getLogger(__name__)
HEADERS = (
    "date",
    "etag",
    "last-modified",
    "content-type",
    "retry-after",
    "x-uniprot-release",
    "x-uniprot-release-date",
)


def timestamp(now):
    from datetime import timezone

    return (
        now()
        .astimezone(timezone.utc)
        .isoformat(timespec="microseconds")
        .replace("+00:00", "Z")
    )


@dataclass
class Audit:
    exchanges: list[HttpExchange] = field(default_factory=list)
    request_count: int = 0
    last_start: float | None = None


@dataclass
class HttpResult:
    data: object
    body: bytes | None
    status: int | None
    headers: dict
    retrieved_at: str
    error: SourceError | None


def request_json(
    url, *, params, role, chemical_id, audit, client, sleep, now, monotonic
):
    audit.request_count += 1
    request_id = audit.request_count
    headers = {
        "Accept": "application/json",
        "User-Agent": f"ligandscout/{version('ligandscout')}",
    }
    final = None
    for attempt in range(1, 4):
        destination, parameters = url, params
        retryable = False
        for hop in range(6):
            delay = (
                max(0, audit.last_start + 0.5 - monotonic())
                if audit.last_start is not None
                else 0
            )
            if delay:
                sleep(delay)
            audit.last_start = monotonic()
            response, transport = None, None
            try:
                response = client.get(
                    destination,
                    params=parameters,
                    headers=headers,
                    timeout=30.0,
                    follow_redirects=False,
                )
            except httpx.RequestError as exc:
                transport = type(exc).__name__
                retryable = isinstance(
                    exc, (httpx.TimeoutException, httpx.NetworkError)
                )
            retrieved = timestamp(now)
            selected = (
                {k: response.headers[k] for k in HEADERS if k in response.headers}
                if response is not None
                else {}
            )
            body = response.content if response is not None else None
            status = response.status_code if response is not None else None
            audit.exchanges.append(
                HttpExchange(
                    sequence=len(audit.exchanges) + 1,
                    request_id=request_id,
                    role=role,
                    chemical_id=chemical_id,
                    attempt=attempt,
                    redirect_hop=hop,
                    url=str(response.request.url.copy_with(query=None))
                    if response is not None
                    else str(destination),
                    parameters=dict(response.request.url.params)
                    if response is not None
                    else dict(parameters)
                    if parameters is not None
                    else dict(httpx.URL(destination).params),
                    retrieved_at=retrieved,
                    http_status=status,
                    transport_error=transport,
                    response_headers=selected,
                    body=body,
                )
            )
            final = HttpResult(None, body, status, selected, retrieved, None)
            if transport:
                final.error = SourceError(
                    "transport", f"{role} request failed ({transport}).", None
                )
                break
            if status in (301, 302, 303, 307, 308):
                location = response.headers.get("location")
                try:
                    target = response.url.join(location) if location else None
                except (httpx.InvalidURL, ValueError):
                    target = None
                if (
                    hop == 5
                    or target is None
                    or target.scheme != "https"
                    or not target.host
                    or target.userinfo
                ):
                    final.error = SourceError(
                        "http",
                        f"{role}: invalid HTTPS redirect or redirect limit.",
                        status,
                    )
                    return final
                destination, parameters = target, None
                continue
            if status == 204 and role == "ccd_search":
                return final
            if status != 200:
                final.error = SourceError(
                    "http", f"{role} returned HTTP {status}.", status
                )
                retryable = status == 429 or 500 <= status <= 599
                break
            try:
                data = strict_json(body) if body.strip() else ""
            except (ValueError, UnicodeError, RecursionError):
                final.error = SourceError(
                    "schema", f"{role} returned invalid JSON.", status
                )
                return final
            if data == "":
                final.error = SourceError(
                    "schema",
                    f"{role} returned an empty response; results are indeterminate.",
                    status,
                )
                retryable = True
                break
            final.data = data
            return final
        if not retryable or attempt == 3:
            return final
        delay = (
            retry_delay(response, attempt - 1, now)
            if response is not None
            else 2 ** (attempt - 1)
        )
        if delay > 60:
            final.error.message += " Retry later: requested wait exceeds 60 seconds."
            return final
        LOG.warning("%s: retry %s/3 in %s seconds", role, attempt + 1, delay)
        sleep(delay)
    raise AssertionError("Unreachable retry state")
