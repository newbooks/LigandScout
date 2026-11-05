"""Stage BindingDB evidence and current-run reports before publication."""

import hashlib
import json
import os
import re
import tempfile
from dataclasses import asdict, replace
from pathlib import Path

from ligandscout.bindingdb import strict_json
from ligandscout.models import as_run, upgrade_result
from ligandscout.report import (
    render_html,
    render_index,
    render_report,
    render_tsv,
    source_slug,
)


def json_bytes(value):
    return (json.dumps(value, indent=2, ensure_ascii=False) + "\n").encode("utf-8")


def raw_extension(body):
    try:
        strict_json(body)
        return "json"
    except (ValueError, UnicodeError, RecursionError):
        return "txt"


def save_results(outcome, *, root: Path) -> Path:
    run = as_run(outcome)
    run = replace(
        run,
        sources=[replace(s, result=upgrade_result(s.result)) for s in run.sources],
    )
    folder = root / run.accession
    sources = folder / "sources"
    raw_files, html_files, source_files, manifests, stale = [], [], [], [], []
    for source in run.sources:
        result = source.result
        if result.accession != run.accession:
            raise ValueError("Cannot combine different accessions.")
        slug = source_slug(result)
        if result.source == "UniProt":
            entries, referenced = [], set()
            for exchange in source.exchanges:
                entry = asdict(exchange)
                body = entry.pop("body")
                entry["raw_path"], entry["sha256"] = None, None
                if body is not None:
                    relative = (
                        f"uniprot-raw/{exchange.sequence:06d}.{raw_extension(body)}"
                    )
                    entry["raw_path"], entry["sha256"] = (
                        relative,
                        hashlib.sha256(body).hexdigest(),
                    )
                    path = sources / relative
                    referenced.add(path)
                    raw_files.append((path, body))
                entries.append(entry)
            manifests.append(
                (
                    sources / "uniprot-requests.json",
                    json_bytes(
                        {
                            "schema_version": 1,
                            "accession": run.accession,
                            "exchanges": entries,
                        }
                    ),
                )
            )
            raw_dir = sources / "uniprot-raw"
            if raw_dir.is_dir():
                stale.extend(
                    path
                    for path in raw_dir.iterdir()
                    if re.fullmatch(r"[0-9]{6,}\.(json|txt)", path.name)
                    and path.is_file()
                    and path not in referenced
                )
        else:
            raw_path = None
            if source.raw_body is not None:
                raw_path = sources / f"{slug}.raw.{raw_extension(source.raw_body)}"
                raw_files.append((raw_path, source.raw_body))
            stale.extend(
                p
                for p in (sources / f"{slug}.raw.json", sources / f"{slug}.raw.txt")
                if p != raw_path
            )
        html_files.append(
            (folder / f"{slug}-ligands.html", render_html(result).encode("utf-8"))
        )
        source_files.append((sources / f"{slug}.json", json_bytes(asdict(result))))
    files = (
        raw_files
        + html_files
        + [
            (folder / "ligands.tsv", render_tsv(run).encode("utf-8")),
            (folder / "report.md", render_report(run).encode("utf-8")),
            (folder / "index.html", render_index(run).encode("utf-8")),
        ]
        + manifests
        + source_files
    )
    staged = []
    try:
        for path, body in files:
            path.parent.mkdir(parents=True, exist_ok=True)
            with tempfile.NamedTemporaryFile(
                dir=path.parent, prefix=".ligandscout-", delete=False
            ) as handle:
                temp = Path(handle.name)
                staged.append((temp, path))
                handle.write(body)
        for temp, path in staged:
            os.replace(temp, path)
        for path in stale:
            path.unlink(missing_ok=True)
    finally:
        for temp, _ in staged:
            temp.unlink(missing_ok=True)
    return folder
