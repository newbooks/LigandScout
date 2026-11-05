"""Source-independent records and provenance."""

from dataclasses import dataclass


@dataclass
class Measurement:
    compound_id: str
    compound_name: str | None
    smiles: str | None
    inchikey: str | None
    target_description: str | None
    evidence_type: str
    measurement_type: str
    measurement_value: str
    measurement_unit: str | None
    pmid: str | None
    doi: str | None
    source_url: str


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


@dataclass
class SearchOutcome:
    result: SourceResult
    raw_body: bytes | None
