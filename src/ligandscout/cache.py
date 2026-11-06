"""Reuse complete source results only when their saved evidence is intact."""

import hashlib
import logging
from pathlib import Path

from ligandscout.bindingdb import parse_records, strict_json
from ligandscout.models import SearchOutcome, source_result_from_dict

LOG = logging.getLogger(__name__)
SLUGS = {"BindingDB": "bindingdb", "UniProt": "uniprot", "RCSB PDB": "rcsb"}


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _manifest_sources(folder, slug, accession):
    manifest = strict_json((folder / "sources" / f"{slug}-requests.json").read_bytes())
    _require(isinstance(manifest, dict), "invalid request manifest")
    _require(
        type(manifest.get("schema_version")) is int and manifest["schema_version"] == 1,
        "unsupported request manifest",
    )
    _require(manifest.get("accession") == accession, "manifest accession mismatch")
    entries = manifest.get("exchanges")
    _require(isinstance(entries, list) and bool(entries), "missing request exchanges")
    request_ids = set()
    paths = set()
    for sequence, entry in enumerate(entries, 1):
        _require(isinstance(entry, dict), "invalid exchange")
        _require(
            type(entry.get("sequence")) is int and entry["sequence"] == sequence,
            "invalid exchange sequence",
        )
        request_id = entry.get("request_id")
        _require(type(request_id) is int and request_id > 0, "invalid request ID")
        request_ids.add(request_id)
        relative, digest = entry.get("raw_path"), entry.get("sha256")
        status = entry.get("http_status")
        if status is None:
            _require(
                relative is None
                and digest is None
                and isinstance(entry.get("transport_error"), str),
                "invalid transport provenance",
            )
            continue
        _require(type(status) is int and 100 <= status <= 599, "invalid HTTP status")
        # Never follow arbitrary paths from a manifest, including parent traversal/symlinks.
        _require(
            relative
            in (f"{slug}-raw/{sequence:06d}.json", f"{slug}-raw/{sequence:06d}.txt"),
            "unexpected raw response path",
        )
        path = folder / "sources" / relative
        _require(
            not path.is_symlink() and not path.parent.is_symlink(),
            "raw response is a symlink",
        )
        _require(
            isinstance(digest, str)
            and hashlib.sha256(path.read_bytes()).hexdigest() == digest,
            "missing or changed raw response",
        )
        _require(relative not in paths, "repeated raw path")
        paths.add(relative)
    return request_ids


def _check_request_refs(value, request_ids):
    if isinstance(value, dict):
        for key, item in value.items():
            if key == "request_ids":
                _require(
                    isinstance(item, list)
                    and all(type(i) is int and i in request_ids for i in item),
                    "result references an absent request",
                )
            else:
                _check_request_refs(item, request_ids)
    elif isinstance(value, list):
        for item in value:
            _check_request_refs(item, request_ids)


def load_completed(root: Path, accession: str, source: str) -> SearchOutcome | None:
    """Read compatible successful/empty results; return None for incomplete evidence."""
    slug = SLUGS[source]
    folder = root / accession
    normalized = folder / "sources" / f"{slug}.json"
    if not normalized.is_file():
        return None
    try:
        data = strict_json(normalized.read_bytes())
        result = source_result_from_dict(data)
        _require(
            result.accession == accession and result.source == source,
            "source identity mismatch",
        )
        _require(
            result.status in ("success", "no_matches") and result.error is None,
            "previous search did not complete",
        )
        _require(result.counts.get("complete") is True, "incomplete result counts")
        report = folder / f"{slug}-ligands.html"
        _require(
            report.is_file() and report.stat().st_size > 0, "missing source report"
        )
        if source == "BindingDB":
            variants = [
                p
                for p in (
                    folder / "sources" / "bindingdb.raw.json",
                    folder / "sources" / "bindingdb.raw.txt",
                )
                if p.is_file()
            ]
            _require(len(variants) == 1, "missing or ambiguous raw response")
            _require(
                parse_records(variants[0].read_bytes()) == result.records,
                "raw response and normalized records disagree",
            )
        else:
            request_ids = _manifest_sources(folder, slug, accession)
            _check_request_refs(data, request_ids)
    except (
        OSError,
        ValueError,
        TypeError,
        KeyError,
        UnicodeError,
        RecursionError,
    ) as exc:
        LOG.info(
            "Saved %s results for %s are not reusable (%s); querying again",
            source,
            accession,
            exc,
        )
        return None
    LOG.info(
        "Skipping %s for %s: reusing completed results retrieved at %s",
        source,
        accession,
        result.retrieved_at,
    )
    return SearchOutcome(result, None, reused=True)
