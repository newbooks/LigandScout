"""Exact-accession ChEMBL searches with auditable assay selection."""

import json
import logging
import re
from datetime import timezone
from decimal import Decimal
from importlib.metadata import version

import httpx

from ligandscout.bindingdb import retry_delay
from ligandscout.models import (
    HttpExchange,
    Measurement,
    SearchOutcome,
    SourceError,
    SourceResult,
)

BASE = "https://www.ebi.ac.uk/chembl/api/data/"
POLICY = "exact-single-protein-v1"
LOG = logging.getLogger(__name__)
ID_FIELDS = {
    "target": "target_chembl_id",
    "activity": "activity_id",
    "assay": "assay_chembl_id",
    "molecule": "molecule_chembl_id",
    "document": "document_chembl_id",
}
TEXT_ACTIVITY = (
    "assay_type",
    "assay_description",
    "molecule_pref_name",
    "assay_variant_accession",
    "assay_variant_mutation",
    "activity_comment",
    "data_validity_comment",
    "canonical_smiles",
)
FAMILY_FIELDS = ("type", "value", "relation", "units", "upper_value", "text_value")


class InvalidData(ValueError):
    """Invalid or incomplete source data; never interpreted as no matches."""


class ServiceFailure(Exception):
    def __init__(self, kind, message, status=None):
        self.error = SourceError(kind, message, status)
        super().__init__(message)


def timestamp(now):
    return (
        now()
        .astimezone(timezone.utc)
        .isoformat(timespec="microseconds")
        .replace("+00:00", "Z")
    )


def decode(body):
    def reject(_):
        raise InvalidData("Nonfinite JSON constant.")

    data = json.loads(body, parse_float=Decimal, parse_constant=reject)
    if not isinstance(data, dict):
        raise InvalidData("Expected a JSON object.")
    return data


def text(value):
    if value is None:
        return None
    if not isinstance(value, str):
        raise InvalidData("Expected optional text.")
    return value if value.strip() else None


def identifier(value):
    if not isinstance(value, str) or not re.fullmatch(r"CHEMBL[0-9]+", value):
        raise InvalidData("Invalid ChEMBL identifier.")
    return value


def digits(value, *, zero=False):
    if value is None:
        return None
    if type(value) is int:
        value = str(value)
    if not isinstance(value, str) or not re.fullmatch(r"[0-9]+", value):
        raise InvalidData("Expected decimal identifier.")
    if not zero and not value.lstrip("0"):
        raise InvalidData("Identifier out of range.")
    return value


def number(value):
    if value is None or isinstance(value, str):
        return text(value)
    if type(value) is int or isinstance(value, Decimal) and value.is_finite():
        return str(value)
    raise InvalidData("Expected finite numeric value or text.")


def optional_object(row, key):
    value = row.get(key)
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise InvalidData("Expected optional object.")
    return value


def validate(resource, row):
    """Validate consumed fields, including rows later excluded by policy."""
    if not isinstance(row, dict):
        raise InvalidData("Expected resource object.")
    if resource == "activity":
        value = row.get("activity_id")
        if type(value) is not int or value <= 0:
            raise InvalidData("Invalid activity identifier.")
        for key in ("assay_chembl_id", "molecule_chembl_id", "target_chembl_id"):
            identifier(row.get(key))
        for key in TEXT_ACTIVITY:
            text(row.get(key))
        for prefix in ("", "standard_"):
            for key in FAMILY_FIELDS:
                (number if key in ("value", "upper_value") else text)(
                    row.get(prefix + key)
                )
        number(row.get("pchembl_value"))
        duplicate = row.get("potential_duplicate")
        if duplicate is not None and not (
            type(duplicate) is bool or type(duplicate) is int and duplicate in (0, 1)
        ):
            raise InvalidData("Invalid duplicate flag.")
    else:
        identifier(row.get(ID_FIELDS[resource]))
    if row.get("document_chembl_id") is not None:
        identifier(row["document_chembl_id"])
    if resource == "target":
        if not text(row.get("target_type")):
            raise InvalidData("Missing target type.")
        components = row.get("target_components")
        if not isinstance(components, list):
            raise InvalidData("Missing target components.")
        for component in components:
            if not isinstance(component, dict):
                raise InvalidData("Invalid target component.")
            text(component.get("accession"))
            text(component.get("component_type"))
        text(row.get("pref_name"))
        text(row.get("organism"))
        digits(row.get("tax_id"))
    elif resource == "assay":
        identifier(row.get("target_chembl_id"))
        for key in ("assay_type", "description", "relationship_type"):
            text(row.get(key))
        confidence = row.get("confidence_score")
        if confidence is not None and (
            type(confidence) is not int or not 0 <= confidence <= 9
        ):
            raise InvalidData("Invalid assay confidence.")
        variant = optional_object(row, "variant_sequence")
        text(variant.get("accession"))
        text(variant.get("mutation"))
        digits(variant.get("isoform"), zero=True)
    elif resource == "molecule":
        text(row.get("pref_name"))
        structures = optional_object(row, "molecule_structures")
        for key in ("canonical_smiles", "standard_inchi_key"):
            text(structures.get(key))
        parent = optional_object(row, "molecule_hierarchy").get("parent_chembl_id")
        if parent is not None:
            identifier(parent)
    elif resource == "document":
        digits(row.get("pubmed_id"))
        text(row.get("doi"))


class ChemblClient:
    def __init__(self, client, sleep, now, monotonic):
        self.client, self.sleep, self.now, self.monotonic = (
            client,
            sleep,
            now,
            monotonic,
        )
        self.exchanges = []
        self.last_start = None
        self.request_id = 0
        self.last_status = None
        self.observed = {}

    def safe_url(self, url, original, *, paging=False):
        try:
            destination = original.join(url)
            if (
                destination.scheme != "https"
                or destination.host != "www.ebi.ac.uk"
                or destination.port not in (None, 443)
                or destination.userinfo
                or destination.fragment
            ):
                raise ValueError()
            if destination.path.rstrip("/") != original.path.rstrip("/"):
                raise ValueError()
            pairs = list(destination.params.multi_items())
            if len(pairs) != len(dict(pairs)):
                raise ValueError()
            params, expected = dict(pairs), dict(original.params)
            if paging:
                for key in ("offset", "limit"):
                    value = params.pop(key, None)
                    if not isinstance(value, str) or not re.fullmatch(r"[0-9]+", value):
                        raise ValueError()
                    if int(value) < (1 if key == "limit" else 0):
                        raise ValueError()
                    expected.pop(key, None)
            if params != expected:
                raise ValueError()
            return destination
        except (ValueError, httpx.InvalidURL):
            raise InvalidData("Unsafe or changed ChEMBL request URL.") from None

    def get(self, url, role, requested_ids):
        self.request_id += 1
        self.last_status = None
        headers = {
            "Accept": "application/json",
            "User-Agent": f"ligandscout/{version('ligandscout')}",
        }
        for attempt in range(1, 4):
            current = httpx.URL(url)
            retry = False
            for hop in range(6):
                if self.last_start is not None:
                    delay = max(0, self.last_start + 0.5 - self.monotonic())
                    if delay:
                        self.sleep(delay)
                self.last_start = self.monotonic()
                response = None
                transport = None
                try:
                    response = self.client.get(
                        current, headers=headers, timeout=30.0, follow_redirects=False
                    )
                except httpx.RequestError as exc:
                    transport = type(exc).__name__
                    retry = isinstance(
                        exc, (httpx.TimeoutException, httpx.NetworkError)
                    )
                self.last_status = (
                    response.status_code if response is not None else None
                )
                self.exchanges.append(
                    HttpExchange(
                        len(self.exchanges) + 1,
                        self.request_id,
                        role,
                        sorted(requested_ids),
                        attempt,
                        hop,
                        str(current),
                        dict(current.params),
                        timestamp(self.now),
                        self.last_status,
                        transport,
                        response.headers.get("Location")
                        if response is not None
                        else None,
                        response.content if response is not None else None,
                    )
                )
                if response is None:
                    error = ServiceFailure(
                        "transport", f"ChEMBL request failed ({transport})."
                    )
                    delay = 2 ** (attempt - 1)
                    break
                if response.status_code in (301, 302, 303, 307, 308):
                    try:
                        location = response.headers.get("Location")
                        if not location or hop == 5:
                            raise InvalidData("Redirect limit or missing Location.")
                        current = self.safe_url(location, current)
                    except InvalidData:
                        raise ServiceFailure(
                            "http",
                            "ChEMBL redirect policy failed.",
                            response.status_code,
                        ) from None
                    continue
                if response.status_code == 200:
                    try:
                        return decode(response.content)
                    except (ValueError, UnicodeError, RecursionError):
                        raise ServiceFailure(
                            "schema", "ChEMBL returned invalid JSON or envelope.", 200
                        ) from None
                error = ServiceFailure(
                    "http",
                    f"ChEMBL returned HTTP {response.status_code}.",
                    response.status_code,
                )
                retry = (
                    response.status_code == 429 or 500 <= response.status_code <= 599
                )
                delay = retry_delay(response, attempt - 1, self.now)
                break
            if not retry or attempt == 3:
                raise error
            if delay > 60:
                raise ServiceFailure(
                    "http",
                    "ChEMBL requested wait exceeds 60 seconds; retry later.",
                    self.last_status,
                )
            LOG.warning("ChEMBL: retry %s/3 in %s seconds", attempt + 1, delay)
            self.sleep(delay)
        raise AssertionError("Unreachable retry state")

    def collection(self, resource, params, role, requested_ids):
        params = {
            **params,
            "limit": "100",
            "offset": "0",
            "order_by": ID_FIELDS[resource],
        }
        original = httpx.URL(BASE + resource + ".json", params=params)
        url, seen, records, total = original, set(), {}, None
        while True:
            if str(url) in seen:
                raise InvalidData("Repeated pagination URL.")
            seen.add(str(url))
            data = self.get(url, role, requested_ids)
            rows, meta = (
                data.get("activities" if resource == "activity" else resource + "s"),
                data.get("page_meta"),
            )
            if not isinstance(rows, list) or not isinstance(meta, dict):
                raise InvalidData("Missing collection or pagination metadata.")
            for key, minimum in (("offset", 0), ("limit", 1), ("total_count", 0)):
                if type(meta.get(key)) is not int or meta[key] < minimum:
                    raise InvalidData("Invalid pagination number.")
            if meta["offset"] != int(url.params["offset"]) or len(rows) > meta["limit"]:
                raise InvalidData("Inconsistent page offset or size.")
            if total is None:
                total = meta["total_count"]
            if total != meta["total_count"]:
                raise InvalidData(
                    "ChEMBL totals changed; incomplete results, retry later."
                )
            for row in rows:
                validate(resource, row)
                key = str(row[ID_FIELDS[resource]])
                observed = self.observed.setdefault(resource, {})
                if key in observed and observed[key] != row:
                    raise InvalidData("Conflicting resource ID across queries.")
                observed[key] = row
                if key in records and records[key] != row:
                    raise InvalidData("Conflicting duplicate resource ID.")
                records[key] = row
            if len(records) > total:
                raise InvalidData("More distinct records than reported total.")
            following = meta.get("next")
            if (
                "next" not in meta
                or following is not None
                and (not isinstance(following, str) or not following.strip())
            ):
                raise InvalidData("Missing or invalid next page.")
            if following is None:
                if len(records) != total:
                    raise InvalidData("ChEMBL pagination incomplete; retry later.")
                return records
            if not rows:
                raise InvalidData("Empty nonfinal page.")
            next_url = self.safe_url(following, url, paging=True)
            if int(next_url.params["offset"]) != meta["offset"] + len(rows):
                raise InvalidData("Pagination did not advance correctly.")
            url = next_url

    def batches(self, resource, ids, role):
        ids, output = sorted(set(ids)), {}
        for start in range(0, len(ids), 50):
            batch = ids[start : start + 50]
            rows = self.collection(
                resource, {ID_FIELDS[resource] + "__in": ",".join(batch)}, role, batch
            )
            if set(rows) != set(batch):
                raise InvalidData("Missing or unexpected joined entity.")
            output.update(rows)
        return output


def target_reasons(target, accession):
    reasons = []
    components = target["target_components"]
    if target["target_type"] != "SINGLE PROTEIN":
        reasons.append("target_type")
    if len(components) != 1:
        reasons.append("component_count")
    elif text(components[0].get("component_type")) not in (None, "PROTEIN"):
        reasons.append("component_type")
    if not any(c.get("accession") == accession for c in components):
        reasons.append("component_accession_mismatch")
    if not components or any(
        not text(c.get("accession")) or not text(c.get("component_type"))
        for c in components
    ):
        reasons.append("unverified_target_assignment")
    return sorted(reasons)


def assay_reasons(assay):
    reasons = []
    relation, confidence, kind = (
        text(assay.get("relationship_type")),
        assay.get("confidence_score"),
        text(assay.get("assay_type")),
    )
    if relation is None or confidence is None or kind is None:
        reasons.append("unverified_target_assignment")
    if relation == "H":
        reasons.append("homolog_assignment")
    elif relation is not None and relation != "D":
        reasons.append("non_direct_assignment")
    if confidence is not None and confidence != 9:
        reasons.append("confidence_not_9")
    if kind is not None and kind not in ("B", "F"):
        reasons.append("unsupported_assay_type")
    return sorted(reasons)


def joined_metadata(activity, assay):
    if activity["target_chembl_id"] != assay["target_chembl_id"]:
        raise InvalidData("Activity and assay target disagree.")
    if (
        text(activity.get("assay_type"))
        and text(assay.get("assay_type"))
        and activity["assay_type"] != assay["assay_type"]
    ):
        raise InvalidData("Activity and assay type disagree.")
    # ChEMBL activities may originate in a different document from the assay.
    # Cite the measurement's own document, falling back only when it is absent.
    activity_document = activity.get("document_chembl_id")
    assay_document = assay.get("document_chembl_id")
    variant = optional_object(assay, "variant_sequence")
    output = {}
    for key, field in (
        ("accession", "assay_variant_accession"),
        ("mutation", "assay_variant_mutation"),
    ):
        left, right = text(variant.get(key)), text(activity.get(field))
        if left and right and left != right:
            raise InvalidData("Activity and assay variant disagree.")
        output[key] = left or right
    output["isoform"] = digits(variant.get("isoform"), zero=True)
    return activity_document or assay_document, output


def normalize(activity, assay, target, molecule, document, accession):
    document_id, variant = joined_metadata(activity, assay)
    families = {}
    for name, prefix in (("published", ""), ("standard", "standard_")):
        families[name] = {
            key: (number if key in ("value", "upper_value") else text)(
                activity.get(prefix + key)
            )
            for key in FAMILY_FIELDS
        }
    basis = next(
        (
            key
            for key in ("standard", "published")
            if families[key]["type"] is not None and families[key]["value"] is not None
        ),
        "qualitative",
    )
    if basis == "qualitative":
        selected = dict.fromkeys(FAMILY_FIELDS)
        selected["type"] = families["standard"]["type"] or families["published"]["type"]
        selected["text_value"] = (
            families["standard"]["text_value"]
            or families["published"]["text_value"]
            or text(activity.get("activity_comment"))
        )
    else:
        selected = families[basis]
    structures = optional_object(molecule, "molecule_structures")
    return Measurement(
        compound_id=activity["molecule_chembl_id"],
        compound_name=text(molecule.get("pref_name"))
        or text(activity.get("molecule_pref_name")),
        smiles=text(structures.get("canonical_smiles")),
        inchikey=text(structures.get("standard_inchi_key")),
        target_description=text(target.get("pref_name")),
        evidence_type="reported_binding_assay"
        if assay["assay_type"] == "B"
        else "reported_functional_assay",
        measurement_type=selected["type"],
        measurement_value=selected["value"],
        measurement_unit=selected["units"],
        pmid=digits(document.get("pubmed_id")),
        doi=text(document.get("doi")),
        source_url="https://www.ebi.ac.uk/chembl/web_components/explore/compound/"
        + activity["molecule_chembl_id"],
        source_record_id=str(activity["activity_id"]),
        target_id=target["target_chembl_id"],
        target_accessions=[accession],
        assay_id=assay["assay_chembl_id"],
        document_id=document_id,
        assay_type=assay["assay_type"],
        assay_description=text(assay.get("description")),
        target_organism=text(target.get("organism")),
        target_tax_id=digits(target.get("tax_id")),
        target_relationship=assay["relationship_type"],
        target_confidence=assay["confidence_score"],
        target_match_status="exact_accession_direct_assignment",
        variant_accession=variant["accession"],
        variant_mutation=variant["mutation"],
        variant_isoform=variant["isoform"],
        activity_comment=text(activity.get("activity_comment")),
        data_validity_comment=text(activity.get("data_validity_comment")),
        potential_duplicate=None
        if activity.get("potential_duplicate") is None
        else bool(activity["potential_duplicate"]),
        parent_compound_id=optional_object(molecule, "molecule_hierarchy").get(
            "parent_chembl_id"
        ),
        measurement_relation=selected["relation"],
        measurement_upper_value=selected["upper_value"],
        measurement_text=selected["text_value"],
        measurement_basis=basis,
        source_values={
            **families,
            "pchembl_value": number(activity.get("pchembl_value")),
        },
    )


def fetch_chembl(accession, *, client, sleep, now, monotonic):
    api = ChemblClient(client, sleep, now, monotonic)
    warnings, records, error = [], [], None
    targets, eligible, activities = {}, {}, {}
    audit = {"outcome_reason": None, "targets": [], "activities": []}
    LOG.info("Searching ChEMBL for %s", accession)

    def release(role):
        try:
            data = api.get(BASE + "status.json", role, [])
            if data.get("status") is not None and data["status"] != "UP":
                warnings.append(
                    {
                        "code": "release_service_not_up",
                        "message": "ChEMBL status is not UP.",
                        "role": role,
                    }
                )
            release_version = text(data.get("chembl_db_version"))
            if not release_version:
                raise InvalidData("Missing release version.")
            return release_version
        except (ServiceFailure, InvalidData):
            warnings.append(
                {
                    "code": "release_metadata_unavailable",
                    "message": "ChEMBL release metadata unavailable.",
                    "role": role,
                }
            )
            return None

    before = release("status_before")
    try:
        targets = api.collection(
            "target",
            {"target_components__accession": accession},
            "targets",
            [accession],
        )
        for tid, target in sorted(targets.items()):
            reasons = target_reasons(target, accession)
            audit["targets"].append(
                {"target_id": tid, "included": not reasons, "reasons": reasons}
            )
            if not reasons:
                eligible[tid] = target
        for tid in sorted(eligible):
            rows = api.collection(
                "activity", {"target_chembl_id": tid}, "activities", [tid]
            )
            for aid, activity in rows.items():
                if activity["target_chembl_id"] != tid:
                    raise InvalidData("Activity target differs from requested target.")
                if aid in activities and activities[aid] != activity:
                    raise InvalidData("Conflicting activity ID.")
                activities[aid] = activity
        assays = api.batches(
            "assay", [a["assay_chembl_id"] for a in activities.values()], "assays"
        )
        included, documents = [], {}
        for aid, activity in sorted(activities.items()):
            assay = assays[activity["assay_chembl_id"]]
            doc, _ = joined_metadata(activity, assay)
            reasons = assay_reasons(assay)
            audit["activities"].append(
                {
                    "activity_id": aid,
                    "assay_id": assay["assay_chembl_id"],
                    "target_id": assay["target_chembl_id"],
                    "included": not reasons,
                    "reasons": reasons,
                }
            )
            if not reasons:
                included.append(activity)
                documents[aid] = doc
        molecules = api.batches(
            "molecule", [a["molecule_chembl_id"] for a in included], "molecules"
        )
        docs = api.batches(
            "document", [d for d in documents.values() if d], "documents"
        )
        records = [
            normalize(
                a,
                assays[a["assay_chembl_id"]],
                eligible[a["target_chembl_id"]],
                molecules[a["molecule_chembl_id"]],
                docs.get(documents[str(a["activity_id"])], {}),
                accession,
            )
            for a in included
        ]
        records.sort(
            key=lambda r: (r.compound_id, r.target_id, r.assay_id, r.source_record_id)
        )
        if not records:
            audit["outcome_reason"] = (
                "no_targets"
                if not targets
                else "no_eligible_targets"
                if not eligible
                else "no_activities"
                if not activities
                else "all_activities_excluded"
            )
    except ServiceFailure as exc:
        error = exc.error
    except (InvalidData, UnicodeError, RecursionError) as exc:
        LOG.debug("ChEMBL validation failed", exc_info=True)
        error = SourceError(
            "schema", f"ChEMBL response validation failed: {exc}", api.last_status
        )
    after = release("status_after")
    if before and after and before != after and error is None:
        error = SourceError(
            "schema", "ChEMBL release changed during search; retry later.", 200
        )
    if error:
        records = []
        audit["outcome_reason"] = "incomplete"
    counts = {
        "complete": error is None,
        "candidate_targets": len(api.observed.get("target", {})),
        "eligible_targets": len(eligible),
        "candidate_activities": len(api.observed.get("activity", {})),
        "included_activities": len(records) if not error else None,
        "excluded_activities": len(activities) - len(records) if not error else None,
        "unverified_activities": sum(
            "unverified_target_assignment" in a["reasons"] for a in audit["activities"]
        )
        if not error
        else None,
        "unique_ligands": len({r.compound_id for r in records}) if not error else None,
    }
    warnings.sort(key=lambda w: (0 if w["role"] == "status_before" else 1, w["code"]))
    result = SourceResult(
        schema_version=2,
        accession=accession,
        source="ChEMBL",
        retrieved_at=timestamp(now),
        request={
            "url": BASE + "target.json",
            "parameters": {
                "target_components__accession": accession,
                "limit": "100",
                "offset": "0",
                "order_by": "target_chembl_id",
            },
            "manifest": "chembl-requests.json",
        },
        source_version=before if before and before == after else None,
        status="failed" if error else "success" if records else "no_matches",
        records=records,
        error=error,
        policy=POLICY,
        warnings=warnings,
        counts=counts,
        selection_audit=audit,
    )
    LOG.info(
        "ChEMBL %s: %s; %s compounds, %s activities",
        accession,
        result.status,
        counts["unique_ligands"],
        counts["included_activities"],
    )
    return SearchOutcome(result, None, api.exchanges)
