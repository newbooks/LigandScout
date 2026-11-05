"""Save per-source evidence and overwrite compiled outputs."""

import json
import os
import tempfile
from dataclasses import asdict
from pathlib import Path

from ligandscout.bindingdb import strict_json
from ligandscout.models import SearchOutcome
from ligandscout.report import render_html, render_report, render_tsv


def save_results(outcome: SearchOutcome, *, root: Path) -> Path:
    folder = root / outcome.result.accession
    sources = folder / "sources"
    files = []
    raw_path = None
    if outcome.raw_body is not None:
        try:
            strict_json(outcome.raw_body)
            extension = "json"
        except (ValueError, UnicodeError, RecursionError):
            extension = "txt"
        raw_path = sources / f"bindingdb.raw.{extension}"
        files.append((raw_path, outcome.raw_body))
    files.extend(
        [
            (folder / "ligands.tsv", render_tsv(outcome.result).encode("utf-8")),
            (folder / "report.md", render_report(outcome.result).encode("utf-8")),
            (folder / "bindingdb-ligands.html", render_html(outcome.result).encode("utf-8")),
            (
                sources / "bindingdb.json",
                (
                    json.dumps(asdict(outcome.result), indent=2, ensure_ascii=False)
                    + "\n"
                ).encode("utf-8"),
            ),
        ]
    )
    sources.mkdir(parents=True, exist_ok=True)
    staged = []
    try:
        for path, body in files:
            with tempfile.NamedTemporaryFile(
                dir=path.parent, prefix=".ligandscout-", delete=False
            ) as handle:
                temp = Path(handle.name)
                staged.append((temp, path))
                handle.write(body)
        if raw_path:
            temp, path = staged[0]
            os.replace(temp, path)
        for variant in (sources / "bindingdb.raw.json", sources / "bindingdb.raw.txt"):
            if variant != raw_path:
                variant.unlink(missing_ok=True)
        for temp, path in staged:
            if path != raw_path:
                os.replace(temp, path)
    finally:
        for temp, _ in staged:
            temp.unlink(missing_ok=True)
    return folder
