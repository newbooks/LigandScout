"""Source-independent Markdown and TSV reports."""

import csv
import html
import io
import re
from collections import Counter
from dataclasses import asdict
from pathlib import Path
from urllib.parse import quote

from ligandscout.models import SourceResult

COLUMNS = [
    "accession",
    "source",
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
    "retrieved_at",
]

EVIDENCE_NOTE = (
    "API-associated records, not confirmed binding partners. "
    "Target accession and assay evidence have not been independently verified. "
    "An IC50 measurement alone does not establish direct binding."
)


def render_summary(result: SourceResult, *, folder: Path) -> str:
    """Summarize unique API-associated ligands for the requested accession."""
    lines = [
        f"LigandScout: {result.accession}",
        f"Source: {result.source} | Status: {result.status}",
        f"Unique ligands listed in {result.source} for {result.accession}: "
        f"{len({r.compound_id for r in result.records})}",
        f"Measurements: {len(result.records)}",
    ]
    if result.error:
        lines.append(f"Search failed: {result.error.message}")
    elif result.status == "no_matches":
        lines.extend(
            [
                f"No BindingDB records found for {result.accession}.",
                "This does not establish the absence of binding partners.",
            ]
        )
    if result.records:
        for label, values in (
            (
                "Reported targets",
                Counter(r.target_description or "Unspecified" for r in result.records),
            ),
            ("Measurement types", Counter(r.measurement_type for r in result.records)),
        ):
            ordered = sorted(values.items(), key=lambda item: (-item[1], item[0]))
            lines.append(f"{label} (measurement counts):")
            for name, count in ordered[:5]:
                # Keep remote labels to one bounded terminal line.
                name = "".join(char if char.isprintable() else " " for char in name)
                lines.append(f"  {name[:100]}: {count}")
            if len(ordered) > 5:
                lines.append(
                    f"  ... {len(ordered) - 5} additional labels in saved results"
                )
        lines.append(EVIDENCE_NOTE)
    lines.extend(
        [f"Results: {folder}", "Full measurements: ligands.tsv | Details: report.md | Ligands: bindingdb-ligands.html"]
    )
    return "\n".join(lines) + "\n"


def escape(value):
    text = html.escape(str(value), quote=True).replace("\r", " ").replace("\n", " ")
    return re.sub(r"([\\|`*_{}\[\]()!~])", r"\\\1", text)


def references(record):
    links = []
    if record.pmid:
        links.append(
            f"[PMID {escape(record.pmid)}](https://pubmed.ncbi.nlm.nih.gov/{quote(record.pmid, safe='')}/)"
            if re.fullmatch(r"[0-9]+", record.pmid)
            else escape(record.pmid)
        )
    if record.doi:
        links.append(
            f"[DOI](https://doi.org/{quote(record.doi, safe='/')})"
            if record.doi.startswith("10.")
            else escape(record.doi)
        )
    return "; ".join(links) or "—"


def render_report(result: SourceResult) -> str:
    lines = [
        f"# LigandScout: {result.accession}",
        "",
        "- Mode: simple",
        f"- Source: {escape(result.source)}",
        f"- Status: {result.status}",
        f"- Retrieved at: {result.retrieved_at}",
        f"- Unique compounds: {len({r.compound_id for r in result.records})}",
        f"- Measurements: {len(result.records)}",
        "",
    ]
    count = str(len({r.compound_id for r in result.records})) if not result.error else "Unavailable"
    lines.extend([
        "| Database | Ligands Count | Ligands Report |",
        "| --- | --- | --- |",
        f"| {escape(result.source)} | {count} | [View ligands](bindingdb-ligands.html) |",
        "",
    ])
    if result.error:
        lines.extend([f"Search failed: {escape(result.error.message)}", ""])
    elif result.status == "no_matches":
        lines.extend(
            [
                f"No BindingDB records found for {result.accession}.",
                "This does not establish the absence of binding partners.",
                "",
            ]
        )
    if result.records:
        lines.extend([EVIDENCE_NOTE, ""])
    lines += [
        "| Compound | Reported target | Measurement | Value | Unit | References |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for record in result.records:
        cells = [
            f"[{escape(record.compound_id)}]({record.source_url})",
            escape(record.target_description or "Unspecified"),
            escape(record.measurement_type),
            escape(record.measurement_value),
            escape(record.measurement_unit or "unspecified"),
            references(record),
        ]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines) + "\n"


def render_tsv(result: SourceResult) -> str:
    stream = io.StringIO(newline="")
    writer = csv.writer(
        stream, delimiter="\t", quoting=csv.QUOTE_MINIMAL, lineterminator="\n"
    )
    writer.writerow(COLUMNS)
    for record in result.records:
        values = asdict(record) | {
            "accession": result.accession,
            "source": result.source,
            "retrieved_at": result.retrieved_at,
        }
        writer.writerow([values[column] for column in COLUMNS])
    return stream.getvalue()


def render_html(result: SourceResult) -> str:
    """Render one row per source compound, preserving every measurement in details."""
    esc = html.escape
    grouped = {}
    for record in result.records:
        grouped.setdefault(record.compound_id, []).append(record)
    lines = [
        '<!doctype html><html lang="en"><head><meta charset="utf-8">',
        '<meta name="viewport" content="width=device-width, initial-scale=1">',
        f'<title>{esc(result.accession)} — {esc(result.source)} ligands</title>',
        '<style>body{font:16px system-ui;max-width:1200px;margin:32px auto;padding:0 16px}'
        'table{border-collapse:collapse;width:100%}th,td{text-align:left;padding:10px;'
        'border-bottom:1px solid #ddd;vertical-align:top}th{background:#eee}'
        'td{overflow-wrap:anywhere}summary{cursor:pointer}li{margin:8px 0}</style></head><body>',
        f'<h1>{esc(result.accession)}: {esc(result.source)} ligands</h1>',
        f'<p>Status: {esc(result.status)} | Retrieved: {esc(result.retrieved_at)}</p>',
        f'<p>Unique ligands: {len(grouped)} | Measurements: {len(result.records)}</p>',
        f'<p>{esc(EVIDENCE_NOTE)}</p>',
        '<p>One row per compound ID. Expand measurements for reported targets and assay values. '
        'Use your browser’s Find command to locate a compound or target.</p>',
        '<p><a href="report.md">Summary report</a> | <a href="ligands.tsv">Full TSV</a></p>',
    ]
    if result.error:
        lines.append(f'<p>Search failed: {esc(result.error.message)}</p>')
    elif result.status == 'no_matches':
        lines.append('<p>No records found. This does not establish the absence of binding partners.</p>')
    lines.append('<table><thead><tr><th>Compound ID</th><th>Name</th><th>Reported targets</th>'
                 '<th>Measurements</th></tr></thead><tbody>')
    for compound_id, records in sorted(grouped.items()):
        # Build links from identifiers rather than interpolating remote URLs into HTML.
        url = ('https://www.bindingdb.org/rwd/bind/chemsearch/marvin/MolStructure.jsp?monomerid='
               + quote(compound_id, safe=''))
        names = sorted({r.compound_name for r in records if r.compound_name})
        targets = sorted({r.target_description or 'Unspecified' for r in records})
        lines.append(f'<tr><td><a href="{esc(url, quote=True)}">{esc(compound_id)}</a></td>'
                     f'<td>{esc("; ".join(names) or "Not supplied")}</td>'
                     f'<td>{esc("; ".join(targets))}</td><td><details>'
                     f'<summary>{len(records)} measurements</summary><ul>')
        for r in records:
            value = f'{r.measurement_type}: {r.measurement_value} {r.measurement_unit or "(unit unspecified)"}'
            reference = '; '.join(filter(None, [f'PMID {r.pmid}' if r.pmid else None, r.doi]))
            lines.append(f'<li>{esc(value)} — {esc(r.target_description or "Unspecified")}'
                         f'{"; " + esc(reference) if reference else ""}</li>')
        lines.append('</ul></details></td></tr>')
    lines.append('</tbody></table></body></html>')
    return '\n'.join(lines) + '\n'
