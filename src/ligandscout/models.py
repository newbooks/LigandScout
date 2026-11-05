"""Source-independent records and provenance."""

from dataclasses import dataclass, field, replace


@dataclass
class Measurement:
    compound_id: str
    compound_name: str | None
    smiles: str | None
    inchikey: str | None
    target_description: str | None
    evidence_type: str
    measurement_type: str | None
    measurement_value: str | None
    measurement_unit: str | None
    pmid: str | None
    doi: str | None
    source_url: str
    source_record_id: str | None = None
    target_id: str | None = None
    target_accessions: list[str] = field(default_factory=list)
    assay_id: str | None = None
    document_id: str | None = None
    assay_type: str | None = None
    assay_description: str | None = None
    target_organism: str | None = None
    target_tax_id: str | None = None
    target_relationship: str | None = None
    target_confidence: int | None = None
    target_match_status: str = "unverified"
    variant_accession: str | None = None
    variant_mutation: str | None = None
    variant_isoform: str | None = None
    activity_comment: str | None = None
    data_validity_comment: str | None = None
    potential_duplicate: bool | None = None
    parent_compound_id: str | None = None
    measurement_relation: str | None = None
    measurement_upper_value: str | None = None
    measurement_text: str | None = None
    measurement_basis: str | None = None
    source_values: dict = field(default_factory=dict)


@dataclass
class SourceError:
    kind: str
    message: str
    http_status: int | None


@dataclass
class SourceResult:
    schema_version: int
    accession: str
    source: str
    retrieved_at: str
    request: dict
    source_version: str | None
    status: str
    records: list[Measurement]
    error: SourceError | None
    policy: str = "bindingdb-associated-v1"
    warnings: list[dict] = field(default_factory=list)
    counts: dict = field(default_factory=dict)
    selection_audit: dict = field(default_factory=dict)


@dataclass
class SearchOutcome:
    result: SourceResult
    raw_body: bytes | None
    exchanges: list["HttpExchange"] = field(default_factory=list)


@dataclass
class HttpExchange:
    sequence: int
    request_id: int
    role: str
    requested_ids: list[str]
    attempt: int
    redirect_hop: int
    url: str
    parameters: dict[str, str]
    retrieved_at: str
    http_status: int | None
    transport_error: str | None
    location: str | None
    body: bytes | None


@dataclass
class RunOutcome:
    accession: str
    sources: list[SearchOutcome]


def upgrade_bindingdb(result: SourceResult) -> SourceResult:
    """Add shared metadata without changing the BindingDB adapter contract."""
    complete = result.status != "failed"
    return replace(
        result,
        schema_version=2,
        policy="bindingdb-associated-v1",
        counts={
            "complete": complete,
            "candidate_targets": None,
            "eligible_targets": None,
            "candidate_activities": len(result.records),
            "included_activities": len(result.records) if complete else None,
            "excluded_activities": 0 if complete else None,
            "unverified_activities": None,
            "unique_ligands": len({r.compound_id for r in result.records})
            if complete
            else None,
        },
        selection_audit={
            "outcome_reason": (
                "incomplete"
                if not complete
                else "no_activities"
                if not result.records
                else None
            ),
            "targets": [],
            "activities": [],
        },
    )


def as_run(value: RunOutcome | SearchOutcome | SourceResult) -> RunOutcome:
    """Keep source renderers and callers usable while supporting multiple sources."""
    if isinstance(value, RunOutcome):
        return value
    if isinstance(value, SourceResult):
        value = SearchOutcome(value, None)
    return RunOutcome(value.result.accession, [value])


def source_result_from_dict(data: dict) -> SourceResult:
    """Read saved v1 BindingDB or v2 results; never combine saved search runs."""
    import types
    from typing import get_args, get_origin, get_type_hints

    def check(value, annotation):
        origin = get_origin(annotation)
        if origin is types.UnionType:
            return any(check(value, option) for option in get_args(annotation))
        if origin is list:
            return isinstance(value, list) and all(
                check(v, get_args(annotation)[0]) for v in value
            )
        if origin is dict:
            return isinstance(value, dict) and all(
                check(k, get_args(annotation)[0]) and check(v, get_args(annotation)[1])
                for k, v in value.items()
            )
        if annotation is int:
            return type(value) is int
        if annotation is type(None):
            return value is None
        return isinstance(value, annotation)

    if (
        not isinstance(data, dict)
        or type(data.get("schema_version")) is not int
        or data["schema_version"] not in (1, 2)
    ):
        raise ValueError("Unsupported saved-result schema version.")
    legacy = data["schema_version"] == 1
    if legacy and data.get("source") != "BindingDB":
        raise ValueError("Version 1 supports BindingDB only.")
    hints = get_type_hints(Measurement)
    original_fields = (
        "compound_id",
        "compound_name",
        "smiles",
        "inchikey",
        "target_description",
        "evidence_type",
        "measurement_type",
        "measurement_value",
        "measurement_unit",
        "pmid",
        "doi",
        "source_url",
    )
    if not isinstance(data.get("records"), list):
        raise TypeError("Invalid saved records.")
    records = []
    for row in data["records"]:
        expected = original_fields if legacy else hints
        if not isinstance(row, dict) or not set(expected) <= row.keys():
            raise ValueError("Missing saved measurement fields.")
        for key in expected:
            if not check(row[key], hints[key]):
                raise ValueError(f"Invalid saved measurement field: {key}.")
        if legacy and (
            not isinstance(row["measurement_type"], str)
            or not isinstance(row["measurement_value"], str)
        ):
            raise ValueError("Version 1 requires measurement type and value.")
        records.append(Measurement(**{key: row[key] for key in expected}))
    root_hints = get_type_hints(SourceResult)
    original_root = (
        "schema_version",
        "accession",
        "source",
        "retrieved_at",
        "request",
        "source_version",
        "status",
        "records",
        "error",
    )
    expected_root = original_root if legacy else root_hints
    if not set(expected_root) <= data.keys():
        raise ValueError("Missing saved source fields.")
    for key in expected_root:
        if key in ("records", "error"):
            continue
        if not check(data[key], root_hints[key]):
            raise ValueError(f"Invalid saved source field: {key}.")
    error = data["error"]
    if error is not None:
        if not isinstance(error, dict) or not all(
            k in error and check(error[k], t)
            for k, t in get_type_hints(SourceError).items()
        ):
            raise ValueError("Invalid saved source error.")
        error = SourceError(**{key: error[key] for key in get_type_hints(SourceError)})
    if data["source"] not in ("BindingDB", "ChEMBL") or data["status"] not in (
        "success",
        "no_matches",
        "failed",
    ):
        raise ValueError("Unknown source or status.")
    if (
        (data["status"] == "failed") != (error is not None)
        or data["status"] != "success"
        and records
    ):
        raise ValueError("Inconsistent saved status, records or error.")
    result = SourceResult(
        **{key: data[key] for key in expected_root if key not in ("records", "error")},
        records=records,
        error=error,
    )
    if legacy:
        return upgrade_bindingdb(result)
    expected_counts = {
        "complete",
        "candidate_targets",
        "eligible_targets",
        "candidate_activities",
        "included_activities",
        "excluded_activities",
        "unverified_activities",
        "unique_ligands",
    }
    if (
        set(result.counts) != expected_counts
        or type(result.counts["complete"]) is not bool
    ):
        raise ValueError("Invalid saved counts.")
    if any(
        v is not None and (type(v) is not int or v < 0)
        for k, v in result.counts.items()
        if k != "complete"
    ):
        raise ValueError("Invalid saved count value.")
    if set(result.selection_audit) != {"outcome_reason", "targets", "activities"}:
        raise ValueError("Invalid saved audit.")
    if not check(result.selection_audit["outcome_reason"], str | None):
        raise ValueError("Invalid saved outcome reason.")
    for kind in ("targets", "activities"):
        entries = result.selection_audit[kind]
        keys = {"target_id", "included", "reasons"} | (
            {"activity_id", "assay_id"} if kind == "activities" else set()
        )
        if not isinstance(entries, list):
            raise TypeError("Invalid audit list.")
        for entry in entries:
            if (
                not isinstance(entry, dict)
                or set(entry) != keys
                or type(entry["included"]) is not bool
                or not check(entry["reasons"], list[str])
                or any(
                    not isinstance(entry[k], str)
                    for k in keys - {"included", "reasons"}
                )
            ):
                raise ValueError("Invalid audit entry.")
    if result.counts["complete"] != (result.status != "failed"):
        raise ValueError("Saved completeness disagrees with source status.")
    if (
        result.policy
        != {
            "BindingDB": "bindingdb-associated-v1",
            "ChEMBL": "exact-single-protein-v1",
        }[result.source]
    ):
        raise ValueError("Unknown source policy.")
    if result.status == "failed":
        if any(
            result.counts[k] is not None
            for k in (
                "included_activities",
                "excluded_activities",
                "unverified_activities",
                "unique_ligands",
            )
        ):
            raise ValueError("Failed source must not expose complete counts.")
    else:
        if result.counts["included_activities"] != len(records) or result.counts[
            "unique_ligands"
        ] != len({r.compound_id for r in records}):
            raise ValueError("Saved counts disagree with records.")
        if result.status == "success" and not records:
            raise ValueError("Successful source must have records.")
        if result.source == "ChEMBL":
            if any(result.counts[k] is None for k in expected_counts - {"complete"}):
                raise ValueError("Completed ChEMBL counts must be known.")
            if (
                result.counts["candidate_activities"]
                != result.counts["included_activities"]
                + result.counts["excluded_activities"]
            ):
                raise ValueError("Saved selection counts do not balance.")
    for record in records:
        if result.source == "ChEMBL":
            values = record.source_values
            family_keys = {
                "type",
                "value",
                "relation",
                "units",
                "upper_value",
                "text_value",
            }
            if set(values) != {"published", "standard", "pchembl_value"} or not check(
                values["pchembl_value"], str | None
            ):
                raise ValueError("Invalid saved source values.")
            for family in ("published", "standard"):
                if (
                    not isinstance(values[family], dict)
                    or set(values[family]) != family_keys
                    or not all(check(v, str | None) for v in values[family].values())
                ):
                    raise ValueError("Invalid saved measurement family.")
            if record.measurement_basis not in ("published", "standard", "qualitative"):
                raise ValueError("Invalid saved measurement basis.")
            if (
                record.target_match_status != "exact_accession_direct_assignment"
                or record.target_confidence != 9
                or record.target_relationship != "D"
            ):
                raise ValueError("Invalid saved ChEMBL target assignment.")
    for warning in result.warnings:
        if (
            set(warning) != {"code", "message", "role"}
            or not check(warning["code"], str)
            or not check(warning["message"], str)
            or not check(warning["role"], str | None)
        ):
            raise ValueError("Invalid saved warning.")
    return result
