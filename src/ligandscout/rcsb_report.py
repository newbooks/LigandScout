"""Readable structural evidence without inline JSON."""

from html import escape
from urllib.parse import quote

NOTE = "RCSB-reported residue contacts establish structural proximity, not affinity or physiological binding."


def site_entries(record):
    for site in record.binding_sites:
        model = site["target_model_id"]
        heading = f"{site['entry_id']} ligand {site['ligand_asym_id']} → target {site['target_asym_id']} (model {model if model is not None else 'unspecified'})"
        for residue in site["residues"]:
            pdb = f"{residue['comp_id']} label {residue['label_seq_id']}"
            if residue["auth_seq_id"] is not None:
                pdb += f" (author residue {residue['auth_seq_id']}; author chain unavailable)"
            mapped = (
                f"{residue['uniprot_accession']}:{residue['uniprot_position']}"
                if residue["uniprot_position"]
                else f"UniProt unmapped: {residue['mapping_reason']}"
            )
            yield (
                f"{heading}: {pdb}; {mapped}; {residue['distance']:g} Å; sequence {residue['sequence_check']}",
                residue["sequence_check"] == "mismatch",
            )


def sites_text(record):
    return (text for text, _ in site_entries(record))


def visible_warnings(result):
    from ligandscout.report import visible_warnings as filtered

    return filtered(result)


def render_markdown(result):
    from ligandscout.report import escape as md_escape

    lines = [
        f"## {result.source}",
        "",
        f"Status: {result.status}",
        f"Retrieved: {result.retrieved_at}",
        "",
        NOTE,
        "",
    ]
    if result.error:
        lines += ["Search failed: " + md_escape(result.error.message), ""]
    elif not result.records:
        lines += [
            "No structure-associated ligands reported: "
            + result.selection_audit["outcome_reason"]
            + ".",
            "This does not establish the absence of binding partners.",
            "",
        ]
    for warning in visible_warnings(result):
        lines.append("Warning: " + md_escape(warning["message"]))
    for record in result.records:
        lines += [
            "",
            f"### [{md_escape(record.compound_name)}]({record.source_url})",
            "",
            "SMILES: "
            + md_escape(record.smiles or record.structure_provenance["status"]),
            "",
        ]
        lines.extend("- " + md_escape(text) for text in sites_text(record))
    lines += ["", *(md_escape(x) for x in result.source_context["coverage_limits"])]
    return "\n".join(lines) + "\n"


def render_rcsb_html(result):
    e = lambda value: escape(str(value), quote=True)
    lines = [
        '<!doctype html><html lang="en"><head><meta charset="utf-8">',
        '<meta name="viewport" content="width=device-width, initial-scale=1">',
        f"<title>{e(result.accession)} — RCSB PDB ligands</title>",
        "<style>body{font:16px system-ui;max-width:1200px;margin:32px auto;padding:0 16px}table{border-collapse:collapse;width:100%}td,th{padding:10px;border-bottom:1px solid #ddd;text-align:left;vertical-align:top}code{white-space:pre-wrap;overflow-wrap:anywhere;user-select:all}li{margin:8px 0}.table-wrap{overflow-x:auto}.residue-mismatch{background:#fff3cd;color:#664d03;border-left:4px solid #997404;padding:6px}</style></head><body>",
        '<p><a href="index.html">Summary</a> | <a href="sources/rcsb.json">Saved source JSON</a> | <a href="sources/rcsb-requests.json">Request provenance</a></p>',
        f"<h1>{e(result.accession)}: Structure-associated ligands</h1>",
        f"<p>Status: {e(result.status)} | Retrieved: {e(result.retrieved_at)}</p>",
        f"<p>{NOTE}</p>",
        "<p>Residue contacts highlighted in amber and marked sequence mismatch differ from the mapped UniProt residue. One record per CCD ID; distinct ligand instances and contact sites remain separate below. Index records count grouped ligands.</p>",
    ]
    if result.error:
        lines.append(
            f"<p>Search failed: {e(result.error.message)}. Counts unavailable.</p>"
        )
    elif not result.records:
        lines.append(
            f"<p>No structure-associated ligands: {e(result.selection_audit['outcome_reason'])}. This does not establish the absence of binding partners.</p>"
        )
    else:
        lines.append(
            f"<p>Unique structural ligands: {len(result.records)} | Contact sites: {result.counts['contact_sites']}</p>"
        )
    for warning in visible_warnings(result):
        lines.append(f"<p>Warning: {e(warning['message'])}</p>")
    lines += [
        '<div class="table-wrap"><table><thead><tr><th>Ligand</th><th>SMILES (RCSB CCD)</th><th>Supporting structures</th><th>RCSB-reported residue contacts</th></tr></thead><tbody>'
    ]
    for record in result.records:
        prov = record.structure_provenance
        comp = quote(record.compound_id[4:], safe="")
        structure = (
            f"<code>{e(record.smiles)}</code><p>{e(prov['smiles_field'])}{' (nonstereospecific)' if prov['smiles_field'] == 'SMILES' else ''}</p>"
            if record.smiles
            else f"SMILES unavailable: {e(prov['status'])}"
        )
        structure += f'<p><a href="https://data.rcsb.org/rest/v1/core/chemcomp/{comp}">CCD API source</a></p>'
        entries = sorted({s["entry_id"] for s in record.binding_sites})
        evidence = []
        for entry in entries:
            meta = result.source_context["entries"][entry]
            methods = ", ".join(
                str(x.get("method", ""))
                for x in meta.get("exptl", [])
                if isinstance(x, dict)
            )
            resolutions = ", ".join(
                f"{x:g} Å" for x in meta.get("resolution_combined", [])
            )
            details = "; ".join(x for x in (methods, resolutions) if x)
            evidence.append(
                f'<li><a href="https://www.rcsb.org/structure/{quote(entry, safe="")}">{e(entry)}</a> {e(details)} | <a href="https://data.rcsb.org/rest/v1/core/entry/{quote(entry, safe="")}">Entry API</a></li>'
            )
        contacts = "".join(
            f'<li class="residue-mismatch">{e(text)}</li>'
            if mismatch
            else f"<li>{e(text)}</li>"
            for text, mismatch in site_entries(record)
        )
        links = "".join(
            f'<li><a href="https://data.rcsb.org/rest/v1/core/nonpolymer_entity_instance/{quote(s["entry_id"], safe="")}/{quote(s["ligand_asym_id"], safe="")}">Ligand instance {e(s["entry_id"])}.{e(s["ligand_asym_id"])}</a></li>'
            for s in record.binding_sites
        )
        lines.append(
            f'<tr><td><a href="{e(record.source_url)}">{e(record.compound_name)}</a><p><a href="https://www.rcsb.org/ligand/{comp}">{e(record.compound_id)}</a></p></td><td>{structure}</td><td><ul>{"".join(evidence)}</ul></td><td><details open><summary>Contact sites</summary><ul class="binding-sites">{contacts}</ul><ul>{links}</ul></details></td></tr>'
        )
    lines += [
        "</tbody></table></div>",
        "<h2>Coverage limits</h2><ul>",
        *(f"<li>{e(text)}</li>" for text in result.source_context["coverage_limits"]),
        "</ul></body></html>",
    ]
    return "\n".join(lines) + "\n"
