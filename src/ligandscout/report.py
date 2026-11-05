"""Source-independent Markdown, local HTML and TSV reports."""

import csv
import html
import io
import json
import re
from collections import Counter
from dataclasses import asdict, fields
from pathlib import Path
from urllib.parse import quote

from ligandscout.models import Measurement, SourceResult, as_run

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
] + [f.name for f in fields(Measurement)][12:]

EVIDENCE_NOTE = (
    "API-associated records, not confirmed binding partners. "
    "Target accession and assay evidence have not been independently verified. "
    "An IC50 measurement alone does not establish direct binding."
)
CHEMBL_NOTE = (
    "Exact accession and direct assay-target assignment; physical binding not independently verified. "
    "Counts include inactive, inconclusive, qualitative, and functional assay records."
)


def evidence_note(result):
    return CHEMBL_NOTE if result.source == "ChEMBL" else EVIDENCE_NOTE


def source_slug(result):
    return {"BindingDB": "bindingdb", "ChEMBL": "chembl"}[result.source]


def ligand_count(result):
    return (
        "Unavailable"
        if result.error
        else str(len({r.compound_id for r in result.records}))
    )


def compound_url(result, record):
    identifier = quote(record.compound_id, safe="")
    if result.source == "ChEMBL":
        return (
            "https://www.ebi.ac.uk/chembl/web_components/explore/compound/" + identifier
        )
    return (
        "https://www.bindingdb.org/rwd/bind/chemsearch/marvin/MolStructure.jsp?monomerid="
        + identifier
    )


def terminal_label(value):
    return "".join(char if char.isprintable() else " " for char in str(value))[:100]


def render_summary(value, *, folder: Path) -> str:
    run = as_run(value)
    lines = [f"LigandScout: {run.accession}"]
    for outcome in run.sources:
        result = outcome.result
        lines.extend(
            [
                f"Source: {result.source} | Status: {result.status}",
                f"Unique ligands listed in {result.source} for {result.accession}: {ligand_count(result)}",
                f"Measurements: {'Unavailable' if result.error else len(result.records)}",
            ]
        )
        if result.error:
            lines.append(f"Search failed: {result.error.message}")
        elif result.status == "no_matches":
            if result.source == "BindingDB":
                lines.append(f"No BindingDB records found for {result.accession}.")
            else:
                lines.append(
                    f"No eligible ChEMBL records for {result.accession}: {result.selection_audit.get('outcome_reason')}."
                )
            lines.append("This does not establish the absence of binding partners.")
        if result.records:
            for label, counts in (
                (
                    "Reported targets",
                    Counter(
                        r.target_description or "Unspecified" for r in result.records
                    ),
                ),
                (
                    "Measurement types",
                    Counter(
                        r.measurement_type or "Unspecified" for r in result.records
                    ),
                ),
            ):
                ordered = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
                lines.append(f"{label} (measurement counts):")
                lines.extend(
                    f"  {terminal_label(name)}: {count}" for name, count in ordered[:5]
                )
                if len(ordered) > 5:
                    lines.append(
                        f"  ... {len(ordered) - 5} additional labels in saved results"
                    )
        if result.source == "ChEMBL":
            lines.append(
                f"Policy: {result.policy} | Release: {result.source_version or 'unknown'}"
            )
            lines.append(
                f"Excluded activities: {result.counts.get('excluded_activities') if not result.error else 'Unavailable'}; "
                f"unverified: {result.counts.get('unverified_activities') if not result.error else 'Unavailable'}"
            )
        if result.records or result.source == "ChEMBL":
            lines.append(evidence_note(result))
        lines.extend(
            f"Warning: {terminal_label(w['message'])}" for w in result.warnings
        )
        lines.append("")
    lines.extend(
        [
            f"Results: {folder}",
            "Full measurements: ligands.tsv | Details: report.md | Ligands: "
            + ", ".join(source_slug(s.result) + "-ligands.html" for s in run.sources),
        ]
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


def measurement_value(record):
    if (
        record.measurement_basis is None
    ):  # BindingDB already preserves qualifiers in its value.
        return record.measurement_value or "unspecified"
    value = (
        record.measurement_value
        if record.measurement_value is not None
        else "unspecified"
    )
    if record.measurement_relation:
        value = record.measurement_relation + " " + value
    if record.measurement_upper_value is not None:
        value += " to " + record.measurement_upper_value
    if record.measurement_text:
        value += "; " + record.measurement_text
    return value


def render_report(value) -> str:
    run = as_run(value)
    lines = [
        f"# LigandScout: {run.accession}",
        "",
        "- Mode: simple",
        "",
        "| Database | Ligands Count | Ligands Report |",
        "| --- | --- | --- |",
    ]
    for outcome in run.sources:
        r = outcome.result
        lines.append(
            f"| {r.source} | {ligand_count(r)} | [View ligands]({source_slug(r)}-ligands.html) |"
        )
    lines.append("")
    if len(run.sources) > 1:
        lines.extend(
            [
                (
                    "ChEMBL uses stricter assay-target selection than BindingDB; source counts are not directly equivalent. "
                    "Compounds and measurements can overlap across databases; counts are not a globally unique total."
                ),
                "",
            ]
        )
    for outcome in run.sources:
        result = outcome.result
        lines.extend(
            [
                f"## {result.source}",
                "",
                f"- Source: {result.source}",
                f"- Status: {result.status}",
                f"- Retrieved at: {result.retrieved_at}",
                f"- Unique compounds: {ligand_count(result)}",
                f"- Measurements: {'Unavailable' if result.error else len(result.records)}",
                "",
            ]
        )
        if result.error:
            lines.extend([f"Search failed: {escape(result.error.message)}", ""])
        elif result.status == "no_matches":
            lines.extend(
                [
                    f"No {result.source} records found for {result.accession}.",
                    "This does not establish the absence of binding partners.",
                    "",
                ]
            )
        if result.source == "ChEMBL":
            lines.extend(
                [
                    f"- Policy: {result.policy}",
                    f"- Release: {escape(result.source_version or 'unknown')}",
                    f"- Outcome reason: {result.selection_audit.get('outcome_reason') or 'none'}",
                    f"- Excluded activities: {result.counts.get('excluded_activities') if not result.error else 'Unavailable'}",
                    f"- Unverified activities: {result.counts.get('unverified_activities') if not result.error else 'Unavailable'}",
                    "- [Full selection audit](sources/chembl.json)",
                    "",
                ]
            )
        lines.extend([evidence_note(result), ""])
        lines.extend(f"Warning: {escape(w['message'])}" for w in result.warnings)
        lines.extend(
            [
                "",
                "| Compound | Reported target | Measurement | Value | Unit | References |",
                "| --- | --- | --- | --- | --- | --- |",
            ]
        )
        for record in result.records:
            cells = [
                f"[{escape(record.compound_id)}]({compound_url(result, record)})",
                escape(record.target_description or "Unspecified"),
                escape(record.measurement_type or "Unspecified"),
                escape(measurement_value(record)),
                escape(record.measurement_unit or "unspecified"),
                references(record),
            ]
            lines.append("| " + " | ".join(cells) + " |")
        lines.append("")
    return "\n".join(lines) + "\n"


def render_tsv(value) -> str:
    run = as_run(value)
    stream = io.StringIO(newline="")
    writer = csv.writer(
        stream, delimiter="\t", quoting=csv.QUOTE_MINIMAL, lineterminator="\n"
    )
    writer.writerow(COLUMNS)
    for outcome in run.sources:
        result = outcome.result
        for record in result.records:
            values = asdict(record) | {
                "accession": result.accession,
                "source": result.source,
                "retrieved_at": result.retrieved_at,
            }
            for key, cell in values.items():
                if isinstance(cell, (list, dict)):
                    values[key] = json.dumps(
                        cell, ensure_ascii=False, sort_keys=True, separators=(",", ":")
                    )
                elif isinstance(cell, bool):
                    values[key] = str(cell).lower()
            writer.writerow([values[column] for column in COLUMNS])
    return stream.getvalue()


def html_references(record):
    links = []
    for label, value, base, valid in (
        (
            "PMID",
            record.pmid,
            "https://pubmed.ncbi.nlm.nih.gov/",
            bool(record.pmid and re.fullmatch(r"[0-9]+", record.pmid)),
        ),
        (
            "DOI",
            record.doi,
            "https://doi.org/",
            bool(record.doi and record.doi.startswith("10.")),
        ),
    ):
        if value:
            content = html.escape(label + " " + value)
            links.append(
                f'<a href="{base}{quote(value, safe="/")}">{content}</a>'
                if valid
                else content
            )
    return "; ".join(links)


def render_html(result: SourceResult) -> str:
    """One row per source compound; expandable details preserve all activities."""
    esc = html.escape
    grouped = {}
    for record in result.records:
        grouped.setdefault(record.compound_id, []).append(record)
    measurements = "Unavailable" if result.error else str(len(result.records))
    lines = [
        '<!doctype html><html lang="en"><head><meta charset="utf-8">',
        '<meta name="viewport" content="width=device-width, initial-scale=1">',
        f"<title>{esc(result.accession)} — {esc(result.source)} ligands</title>",
        (
            "<style>body{font:16px system-ui;max-width:1200px;margin:32px auto;padding:0 16px}"
            "table{border-collapse:collapse;width:100%}th,td{text-align:left;padding:10px;"
            "border-bottom:1px solid #ddd;vertical-align:top}th{background:#eee}"
            "td{overflow-wrap:anywhere}summary{cursor:pointer}li{margin:8px 0}"
            ".smiles{display:block;max-width:32rem;white-space:pre-wrap;"
            "overflow-wrap:anywhere;user-select:all;margin-bottom:8px}</style></head><body>"
        ),
        f"<h1>{esc(result.accession)}: {esc(result.source)} ligands</h1>",
        f"<p>Status: {esc(result.status)} | Retrieved: {esc(result.retrieved_at)}</p>",
        f"<p>Unique ligands: {ligand_count(result)} | Measurements: {measurements}</p>",
        f"<p>{esc(evidence_note(result))}</p>",
        (
            "<p>One row per compound ID. Expand measurements for reported targets and assay values. "
            "Select SMILES text to copy it. Use your browser’s Find command to locate a compound or target.</p>"
        ),
        '<p><a href="report.md">Summary report</a> | <a href="ligands.tsv">Full TSV</a></p>',
    ]
    if result.error:
        lines.append(f"<p>Search failed: {esc(result.error.message)}</p>")
    elif result.status == "no_matches":
        lines.append(
            "<p>No records found. This does not establish the absence of binding partners.</p>"
        )
    if result.source == "ChEMBL":
        lines.extend(
            [
                (
                    f"<p>Policy: {esc(result.policy)} | Release: {esc(result.source_version or 'unknown')} | "
                    f"Outcome: {esc(result.selection_audit.get('outcome_reason') or 'success')}</p>"
                ),
                "<p>Excluded activities: "
                + str(
                    result.counts.get("excluded_activities")
                    if not result.error
                    else "Unavailable"
                )
                + "; unverified: "
                + str(
                    result.counts.get("unverified_activities")
                    if not result.error
                    else "Unavailable"
                )
                + '. <a href="sources/chembl.json">Full selection audit</a></p>',
            ]
        )
    lines.extend(f"<p>Warning: {esc(w['message'])}</p>" for w in result.warnings)
    lines.append(
        "<table><thead><tr><th>Compound ID</th><th>Name</th><th>SMILES</th><th>Reported targets</th>"
        "<th>Measurements</th></tr></thead><tbody>"
    )
    for compound_id, records in sorted(grouped.items()):
        names = sorted({r.compound_name for r in records if r.compound_name})
        targets = sorted({r.target_description or "Unspecified" for r in records})
        smiles = sorted({r.smiles for r in records if r.smiles and r.smiles.strip()})
        smiles_cell = (
            "".join(f'<code class="smiles">{esc(value)}</code>' for value in smiles)
            or "Not supplied"
        )
        lines.append(
            f'<tr><td><a href="{esc(compound_url(result, records[0]), quote=True)}">{esc(compound_id)}</a></td>'
            f"<td>{esc('; '.join(names) or 'Not supplied')}</td>"
            f"<td>{smiles_cell}</td>"
            f"<td>{esc('; '.join(targets))}</td><td><details>"
            f"<summary>{len(records)} measurements</summary><ul>"
        )
        for r in records:
            value = f"{r.measurement_type or 'Unspecified'}: {measurement_value(r)} {r.measurement_unit or '(unit unspecified)'}"
            lines.append(
                f"<li>{esc(value)} — {esc(r.target_description or 'Unspecified')} {html_references(r)}"
            )
            metadata = {
                "Activity ID": r.source_record_id,
                "Assay type": r.assay_type,
                "Assay": r.assay_description,
                "Organism": r.target_organism,
                "Taxonomy ID": r.target_tax_id,
                "Target assignment": r.target_match_status,
                "Relationship": r.target_relationship,
                "Confidence": r.target_confidence,
                "Accessions": "; ".join(r.target_accessions),
                "Variant accession": r.variant_accession,
                "Mutation": r.variant_mutation,
                "Isoform": r.variant_isoform,
                "Activity comment": r.activity_comment,
                "Data validity": r.data_validity_comment,
                "Potential duplicate": r.potential_duplicate,
                "Parent compound ID": r.parent_compound_id,
            }
            for key, item in metadata.items():
                if item is not None and item != "":
                    lines.append(f"<div>{esc(key)}: {esc(str(item))}</div>")
            for kind, identifier in (
                ("target", r.target_id),
                ("assay", r.assay_id),
                ("document", r.document_id),
            ):
                if (
                    result.source == "ChEMBL"
                    and identifier
                    and re.fullmatch(r"CHEMBL[0-9]+", identifier)
                ):
                    lines.append(
                        f'<div><a href="https://www.ebi.ac.uk/chembl/web_components/explore/{kind}/{identifier}">{kind}: {identifier}</a></div>'
                    )
            lines.append("</li>")
        lines.append("</ul></details></td></tr>")
    lines.append("</tbody></table></body></html>")
    return "\n".join(lines) + "\n"
