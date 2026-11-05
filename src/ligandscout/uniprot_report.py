"""Source-specific annotation reports; no synthetic affinity measurements."""

import html
import json
import re
from urllib.parse import quote

from ligandscout.report import escape, evidence_note


def site_label(site):
    def position(value):
        number = str(value["value"]) if value["value"] is not None else "?"
        return (
            number
            if value["modifier"] == "EXACT"
            else f"{number} ({value['modifier']})"
        )

    start, end = site["location"]["start"], site["location"]["end"]
    label = position(start) if start == end else position(start) + "–" + position(end)
    if site["residue_sequence"]:
        label += " (" + site["residue_sequence"] + ")"
    if site["sequence_scope"] == "alternate":
        label += f" — {site['location']['sequence']}: alternate sequence; residue letters unavailable"
    products = site["processed_products"]
    if products:
        label += "; products: " + ", ".join(
            p["name"] or p["annotation_path"] for p in products
        )
    if site["ligand_part"]:
        part = site["ligand_part"]
        label += "; ligand part: " + (part["name"] or part["id"] or "Unspecified")
    return label


def evidence_html(items):
    if not items:
        return "<h4>Evidence</h4><p>Evidence not supplied.</p>"
    lines = ["<h4>Evidence</h4><ul>"]
    for item in items:
        label = " / ".join(
            str(item[k])
            for k in ("category", "evidence_code", "source", "source_id")
            if item[k]
        )
        source_id = item["source_id"] or ""
        url = None
        if item["source"] == "PubMed" and re.fullmatch(r"[0-9]+", source_id):
            url = "https://pubmed.ncbi.nlm.nih.gov/" + source_id + "/"
        elif item["source"] == "UniProtKB" and re.fullmatch(
            r"[A-Z0-9]{6,10}", source_id
        ):
            url = (
                "https://www.uniprot.org/uniprotkb/"
                + quote(source_id, safe="")
                + "/entry"
            )
        content = html.escape(label)
        if url:
            content = f'<a href="{url}">{content}</a>'
        lines.append(f"<li>{content}</li>")
    return "\n".join(lines) + "</ul>"


def pretty(value):
    return json.dumps(value, indent=2, ensure_ascii=False)


def summarize_text(texts, max_sentences=3):
    """Select leading source sentences without paraphrasing or dropping qualifiers."""
    text = " ".join(
        " ".join(value.split())
        for value in texts
        if isinstance(value, str) and value.strip()
    )
    sentences, start = [], 0
    for boundary in re.finditer(r"[.!?]\s+(?=[A-Z])", text):
        end = boundary.start() + 1
        if (
            text[:end]
            .lower()
            .endswith(("e.g.", "i.e.", "dr.", "mr.", "mrs.", "vs.", "fig.", "no."))
        ):
            continue
        sentences.append(text[start:end])
        start = boundary.end()
    if text[start:]:
        sentences.append(text[start:])
    return " ".join(sentences[:max_sentences]), len(sentences) > max_sentences


def note_texts(annotation):
    note = annotation.get("note")
    return (
        [
            item.get("value", "")
            for item in note.get("texts", [])
            if isinstance(item, dict)
        ]
        if isinstance(note, dict)
        else []
    )


def context_texts(kind, data):
    """Read known text-bearing fields; never serialize the context object for HTML."""
    texts = [
        item.get("value", "")
        for item in data.get("texts", [])
        if isinstance(item, dict)
    ]
    if kind == "CATALYTIC ACTIVITY":
        reaction = data.get("reaction", {})
        if reaction.get("name"):
            texts.append("Annotated reaction: " + reaction["name"] + ".")
        if reaction.get("ecNumber"):
            texts.append("Enzyme classification: " + reaction["ecNumber"] + ".")
    elif kind == "COFACTOR":
        names = [c.get("name") for c in data.get("cofactors", []) if c.get("name")]
        if names:
            texts.append("Annotated cofactors: " + ", ".join(names) + ".")
        texts.extend(note_texts({"note": data.get("note")}))
    elif kind == "BIOPHYSICOCHEMICAL PROPERTIES":
        for field, label in (
            ("phDependence", "pH dependence"),
            ("temperatureDependence", "Temperature dependence"),
            ("redoxPotential", "Redox potential"),
            ("absorption", "Absorption"),
            ("kineticParameters", "Kinetic parameters"),
        ):
            value = data.get(field)
            if not isinstance(value, dict):
                continue
            for item in (
                value.get("texts", [])
                if isinstance(value.get("texts", []), list)
                else []
            ):
                if (
                    isinstance(item, dict)
                    and isinstance(item.get("value"), str)
                    and item["value"]
                ):
                    texts.append(label + ": " + item["value"])
            if (
                field == "absorption"
                and value.get("max") is not None
                and isinstance(value.get("unit"), str)
            ):
                texts.append(
                    "Absorption maximum: "
                    + str(value["max"])
                    + " "
                    + value["unit"]
                    + "."
                )
            if field == "kineticParameters":
                for key, metric, unit in (
                    ("michaelisConstants", "Michaelis constant", None),
                    ("maximumVelocities", "Maximum velocity", None),
                ):
                    for item in (
                        value.get(key, [])
                        if isinstance(value.get(key, []), list)
                        else []
                    ):
                        if (
                            not isinstance(item, dict)
                            or not isinstance(item.get("constant"), (str, int, float))
                            or isinstance(item.get("constant"), bool)
                            or not isinstance(item.get("unit", unit), str)
                        ):
                            continue
                        sentence = (
                            f"{metric}: {item['constant']} {item.get('unit', unit)}"
                        )
                        if isinstance(item.get("substrate"), str) and item["substrate"]:
                            sentence += " for " + item["substrate"]
                        texts.append(sentence + ".")
    return texts


def structure_details_html(provenance):
    esc = html.escape
    labels = {
        "available": "Exact chemical-ID match in RCSB CCD.",
        "no_exact_mapping": "No exact chemical-ID match found in RCSB CCD.",
        "no_descriptor": "The matched component has no usable SMILES descriptor.",
        "ambiguous": "Candidate components have conflicting structures; no SMILES was selected.",
        "unavailable": "The structure lookup could not be completed.",
        "not_applicable": "No chemical identifier is available for an exact structure lookup.",
    }
    lines = [
        "<h4>Structure information</h4><p>" + esc(labels[provenance["status"]]) + "</p>"
    ]
    if provenance.get("requested_chemical_id"):
        lines.append(
            "<p>Requested chemical identity: "
            + esc(provenance["requested_chemical_id"])
            + ".</p>"
        )
    if provenance["candidate_components"]:
        lines.append("<ul>")
        for candidate in provenance["candidate_components"]:
            cid = candidate["ccd_id"]
            label = (
                "exact identifier match"
                if candidate.get("mapping")
                else "identifier does not match"
            )
            lines.append(
                f'<li><a href="https://www.rcsb.org/ligand/{quote(cid, safe="")}">CCD {esc(cid)}</a>: {label}.'
            )
            if candidate.get("smiles") and provenance["status"] != "available":
                lines.append(f'<code class="smiles">{esc(candidate["smiles"])}</code>')
            lines.append("</li>")
        lines.append("</ul>")
    return "\n".join(lines)


def render_uniprot_html(result):
    esc = html.escape
    context = result.source_context
    lines = [
        '<!doctype html><html lang="en"><head><meta charset="utf-8">',
        '<meta name="viewport" content="width=device-width, initial-scale=1">',
        f"<title>{esc(result.accession)} — UniProt annotations</title>",
        "<style>body{font:16px system-ui;max-width:1400px;margin:32px auto;padding:0 16px}table{border-collapse:collapse;width:100%}th,td{text-align:left;padding:10px;border-bottom:1px solid #ddd;vertical-align:top}td{overflow-wrap:anywhere}.smiles{display:block;white-space:pre-wrap;overflow-wrap:anywhere;user-select:all}.table-wrap{overflow-x:auto}</style></head><body>",
        f"<h1>{esc(result.accession)}: UniProt-annotated ligands and cofactors</h1>",
        '<p><a href="index.html">Search summary</a> | <a href="ligands.tsv">Full TSV</a> | <a href="sources/uniprot.json">Complete annotation data (JSON)</a></p>',
        f"<p>Status: {esc(result.status)} | Retrieved: {esc(result.retrieved_at)}</p>",
        f"<p>{esc(evidence_note(result))}</p>",
        "<p>Binding-site coordinates are 1-based, inclusive UniProt sequence positions, not PDB residue numbers.</p>",
    ]
    if result.error:
        lines.append(f"<p>Search failed: {esc(result.error.message)}</p>")
    else:
        entry = context["entry"]
        primary = entry["primary_accession"]
        source_url = (
            f"https://www.uniprot.org/uniprotkb/{quote(primary, safe='')}/entry"
        )
        api_url = f"https://rest.uniprot.org/uniprotkb/{quote(primary, safe='')}.json"
        source_link = f'<a href="{source_url}">UniProt {esc(primary)}</a>'
        lines.extend(
            [
                "<h2>Protein</h2>",
                f"<p>{esc(entry['protein_name'])} — {esc(entry['organism'])} (taxon {entry['taxon_id']})</p>",
                f"<p>Query: {esc(result.accession)} | Resolved primary accession: {esc(primary)} | {esc(entry['entry_type'])} | Sequence length: {entry['sequence_length']}</p>",
                f"<p>Entry version: {esc(str(entry['entry_audit'].get('entryVersion') or 'unknown'))} | Sequence version: {esc(str(entry['entry_audit'].get('sequenceVersion') or 'unknown'))}</p>",
                f'<p>Annotation source: {source_link} | <a href="{api_url}">UniProt API record</a>. The saved annotation data preserves the retrieved version.</p>',
                f"<p>Identified ligand groups: {len(result.records)} | With chemical IDs: {result.counts['identified_chemical_ids']} | With SMILES: {result.counts['smiles_available']}</p>",
            ]
        )
        if not result.records:
            lines.append(
                "<p>No specifically identified ligands/cofactors annotated in UniProt. This does not establish the absence of binding partners.</p>"
            )
        if result.records:
            lines.extend(
                [
                    "<h2>Annotated ligands and cofactors</h2>",
                    f"<p>Ligand and binding-site source: {source_link}.</p>",
                    '<div class="table-wrap"><table><thead><tr><th>Ligand</th><th>Role</th><th>SMILES / structure provider</th><th>Binding sites</th></tr></thead><tbody>',
                ]
            )
        for record in result.records:
            prov = record.structure_provenance
            chemical = esc(record.compound_id)
            if re.fullmatch(r"CHEBI:[1-9][0-9]*", record.compound_id):
                chemical = f'<a href="https://www.ebi.ac.uk/chebi/{quote(record.compound_id, safe=":")}">{chemical}</a>'
            structure = (
                f'<code class="smiles">{esc(record.smiles)}</code><p>RCSB CCD — {esc(prov["smiles_field"])}'
                + (
                    " (non-stereospecific fallback)"
                    if prov["smiles_field"] == "SMILES"
                    else " (stereochemistry as supplied)"
                )
                + "</p>"
                if record.smiles
                else "<p>SMILES unavailable.</p>"
            )
            sites = (
                '<ul class="binding-sites">'
                + "".join(
                    f"<li>{esc(site_label(site))}</li>" for site in record.binding_sites
                )
                + "</ul>"
                if record.binding_sites
                else "Binding site not annotated."
            )
            lines.append(
                f"<tr><td>{esc(record.compound_name or 'Name not supplied')}<br>{chemical}</td>"
                f"<td>{esc(', '.join(record.ligand_roles))}</td>"
                f"<td>{structure}{structure_details_html(prov)}</td><td>{sites}</td></tr>"
            )
        if result.records:
            lines.append("</tbody></table></div>")
        if context["generic_sites"]:
            lines.append("<h2>Generic binding sites (no identified molecule)</h2>")
            lines.append('<ul class="binding-sites">')
            for binding in context["generic_sites"]:
                lines.append(f"<li>{esc(site_label(binding))}</li>")
            lines.append("</ul>")
        if context["generic_cofactors"]:
            lines.append("<h2>Unidentified cofactors</h2>")
            for annotation in context["generic_cofactors"]:
                summary, _ = summarize_text(note_texts(annotation))
                lines.append(
                    f"<p>{esc(summary or 'A cofactor is annotated without a specific chemical identity.')}</p>"
                )
                lines.append(evidence_html(annotation["evidence"]))
        grouped = {}
        for comment in context["context_comments"]:
            grouped.setdefault(comment["comment_type"], []).extend(
                context_texts(comment["comment_type"], comment["data"])
            )
        for kind, texts in grouped.items():
            summary, truncated = summarize_text(texts)
            if not summary:
                continue
            lines.append(f"<h2>{esc(kind.replace('_', ' ').title())} context</h2>")
            lines.append(f"<p>{esc(summary)}</p>")
            if kind == "CATALYTIC ACTIVITY":
                lines.append(
                    "<p>Reaction participants are context, not automatically binding partners.</p>"
                )
            note = (
                "First three source sentences; additional information is available in the full annotation."
                if truncated
                else "Text extracted from the source annotation."
            )
            lines.append(f"<p>{note} Source: {source_link}.</p>")
    if result.warnings:
        lines.append("<h2>Warnings</h2><ul>")
        lines.extend(f"<li>{esc(w['message'])}</li>" for w in result.warnings)
        lines.append("</ul>")
    return "\n".join(lines) + "\n</body></html>\n"


def render_markdown(result):
    lines = [
        "## UniProt",
        "",
        f"- Status: {result.status}",
        f"- Retrieved: {result.retrieved_at}",
        "- [Annotated ligands, binding sites and evidence](uniprot-ligands.html)",
        "",
        evidence_note(result),
        "",
    ]
    if result.error:
        lines.append("Search failed: " + escape(result.error.message))
    else:
        entry = result.source_context["entry"]
        lines.append(
            f"Protein: {escape(entry['protein_name'])}; primary accession: {entry['primary_accession']}; sequence length: {entry['sequence_length']}."
        )
        if not result.records:
            lines.append(
                "No specifically identified ligands/cofactors annotated in UniProt. Context below may still be available."
            )
        for record in result.records:
            lines.extend(
                [
                    "",
                    f"### {escape(record.compound_name or record.compound_id)}",
                    "",
                    f"- Chemical ID: {escape(record.compound_id)}",
                    f"- Roles: {', '.join(record.ligand_roles)}",
                    f"- SMILES: {escape(record.smiles) if record.smiles else 'unavailable: ' + record.structure_provenance['status']}",
                    "- Sites: "
                    + (
                        "; ".join(escape(site_label(s)) for s in record.binding_sites)
                        or "Binding site not annotated."
                    ),
                ]
            )
            details = {
                "binding_sites": record.binding_sites,
                "annotation_evidence": record.annotation_evidence,
                "cofactor_annotations": record.source_values["cofactor_annotations"],
                "structure_provenance": record.structure_provenance,
            }
            lines.extend(["", "<pre>" + html.escape(pretty(details)) + "</pre>"])
        lines.extend(
            [
                "",
                "### Annotation context",
                "",
                "<pre>" + html.escape(pretty(result.source_context)) + "</pre>",
            ]
        )
    lines.extend("Warning: " + escape(w["message"]) for w in result.warnings)
    return "\n".join(lines) + "\n"
