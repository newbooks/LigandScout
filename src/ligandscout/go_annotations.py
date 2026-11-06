"""UniProt GO terms with explicit, audited GOA relations from QuickGO."""

import html
import logging
import re

from ligandscout.api_http import request_json

LOG = logging.getLogger(__name__)
ASPECTS = {"F": "MF", "P": "BP", "C": "CC"}
QUICKGO_ASPECTS = {
    "molecular_function": "MF",
    "biological_process": "BP",
    "cellular_component": "CC",
}
URL = "https://www.ebi.ac.uk/QuickGO/services/annotation/search"


def uniprot_terms(data):
    terms = {}
    for ref in data.get("uniProtKBCrossReferences", []):
        if ref.get("database") != "GO":
            continue
        go_id = ref.get("id", "")
        if not re.fullmatch(r"GO:[0-9]{7}", go_id):
            raise ValueError("Invalid UniProt GO identifier")
        term = next(
            (
                p.get("value")
                for p in ref.get("properties", [])
                if p.get("key") == "GoTerm"
            ),
            None,
        )
        if not isinstance(term, str) or term[:1] not in ASPECTS or term[1:2] != ":":
            raise ValueError("Invalid UniProt GO term/aspect")
        terms[go_id] = {
            "id": go_id,
            "name": term[2:],
            "aspect": ASPECTS[term[0]],
            "relations": [],
        }
    return sorted(terms.values(), key=lambda t: (t["aspect"], t["id"]))


def enrich_go(result, data, **kwargs):
    """Optional enrichment never turns a completed ligand search into a failure."""
    section = {"status": "unavailable", "terms": [], "request_ids": [], "error": None}
    result.source_context["go_annotations"] = section
    try:
        section["terms"] = uniprot_terms(data)
        if not section["terms"]:
            section["status"] = "no_matches"
            return
        terms = {t["id"]: t for t in section["terms"]}
        accession = result.source_context["entry"]["primary_accession"]
        relations = {key: set() for key in terms}
        page, total_hits, seen = 1, None, 0
        while True:
            response = request_json(
                URL,
                params={
                    "geneProductId": f"UniProtKB:{accession}",
                    "limit": 200,
                    "page": page,
                },
                role="go_annotations",
                chemical_id=None,
                **kwargs,
            )
            section["request_ids"].append(kwargs["audit"].request_count)
            if response.error:
                raise ValueError(response.error.message)
            payload = response.data
            if not isinstance(payload, dict):
                raise ValueError("Invalid QuickGO response")
            rows, hits, info = (
                payload.get("results"),
                payload.get("numberOfHits"),
                payload.get("pageInfo"),
            )
            if (
                not isinstance(rows, list)
                or type(hits) is not int
                or hits < 0
                or not isinstance(info, dict)
            ):
                raise ValueError("Invalid QuickGO pagination")
            pages = info.get("total")
            if type(pages) is not int or pages < 0 or info.get("current") != page:
                raise ValueError("Invalid QuickGO page number")
            if total_hits is not None and total_hits != hits:
                raise ValueError("QuickGO annotation count changed during pagination")
            total_hits = hits
            for row in rows:
                if (
                    not isinstance(row, dict)
                    or row.get("geneProductId") != f"UniProtKB:{accession}"
                ):
                    raise ValueError("QuickGO returned a different gene product")
                term = terms.get(row.get("goId"))
                if term is None:
                    continue
                if QUICKGO_ASPECTS.get(row.get("goAspect")) != term["aspect"]:
                    raise ValueError("QuickGO and UniProt GO aspects disagree")
                relation = row.get("qualifier")
                if relation is not None:
                    if not isinstance(relation, str):
                        raise ValueError("Invalid QuickGO relation")
                    if relation.strip():
                        relations[term["id"]].add(relation)
            seen += len(rows)
            LOG.info(
                "UniProt GO %s: retrieved %s/%s annotations", accession, seen, hits
            )
            if page >= pages:
                if seen != hits:
                    raise ValueError("Incomplete QuickGO annotation pages")
                break
            if not rows or seen >= hits:
                raise ValueError("Incomplete QuickGO pagination")
            page += 1
        for key, term in terms.items():
            term["relations"] = sorted(relations[key])
        section["status"] = "success"
    except (ValueError, TypeError, KeyError, AttributeError) as exc:
        section["error"] = str(exc)
        LOG.warning("UniProt GO annotations unavailable: %s", exc)


def render_go(result):
    e = html.escape
    entry_type = (
        result.source_context.get("entry", {}).get("entry_type") if result else None
    )
    if entry_type == "UniProtKB reviewed (Swiss-Prot)":
        label = "UniProt | SwissProt"
    elif entry_type == "UniProtKB unreviewed (TrEMBL)":
        label = 'UniProt | <s title="Unreviewed UniProtKB/TrEMBL entry">SwissProt</s>'
    else:
        label = "UniProt | SwissProt status unavailable"
    lines = [f'<section id="go-annotations"><h2>GO annotations ({label})</h2>']
    section = result.source_context.get("go_annotations") if result else None
    if not section:
        return "\n".join(
            lines
            + [
                "<p>GO annotations unavailable in saved results; refresh UniProt with --force.</p></section>"
            ]
        )
    lines.append(
        f"<p>Retrieved: {e(result.retrieved_at)}. SwissProt status describes the UniProt entry review status, not the evidence level of each GO annotation. Terms from UniProt; explicit relation labels from QuickGO/GOA. Missing relations are not defaulted.</p>"
    )
    if section["status"] == "unavailable":
        lines.append(
            f"<p>GO relation lookup unavailable: {e(section['error'] or 'Unknown error')}. Available UniProt terms are listed below.</p>"
        )
    for aspect, title in (
        ("MF", "Molecular Function"),
        ("BP", "Biological Process"),
        ("CC", "Cellular Component"),
    ):
        lines.append(f"<h3>{aspect} — {title}</h3>")
        lines.append(
            '<div class="table-wrap"><table class="go-table">'
            '<colgroup><col style="width:72%"><col style="width:28%"></colgroup>'
            '<thead><tr><th scope="col">Description</th><th scope="col">Relation</th></tr></thead><tbody>'
        )
        terms = [t for t in section["terms"] if t["aspect"] == aspect]
        for term in terms:
            relation = "; ".join(term["relations"]) or "not supplied"
            lines.append(
                f'<tr><td><a href="https://www.ebi.ac.uk/QuickGO/term/{term["id"]}">{term["id"]}</a> — {e(term["name"])}</td><td>{e(relation)}</td></tr>'
            )
        lines.append("</tbody></table></div>")
        if not terms:
            lines.append(
                "<p>No GO terms listed.</p>"
                if section["status"] != "unavailable"
                else "<p>No GO terms available.</p>"
            )
    return "\n".join(lines + ["</section>"])


def validate_saved_go(section):
    """Validate the additive GO context while allowing older saved results."""
    if not isinstance(section, dict) or set(section) != {
        "status",
        "terms",
        "request_ids",
        "error",
    }:
        raise ValueError("Invalid saved GO context")
    if section["status"] not in {"success", "no_matches", "unavailable"}:
        raise ValueError("Invalid saved GO status")
    if not isinstance(section["request_ids"], list) or any(
        type(i) is not int or i < 1 for i in section["request_ids"]
    ):
        raise ValueError("Invalid saved GO request references")
    if section["error"] is not None and not isinstance(section["error"], str):
        raise ValueError("Invalid saved GO error")
    if not isinstance(section["terms"], list):
        raise ValueError("Invalid saved GO terms")
    ids = set()
    for term in section["terms"]:
        if (
            not isinstance(term, dict)
            or set(term) != {"id", "name", "aspect", "relations"}
            or not isinstance(term["id"], str)
            or not re.fullmatch(r"GO:[0-9]{7}", term["id"])
            or term["id"] in ids
            or term["aspect"] not in {"MF", "BP", "CC"}
            or not isinstance(term["name"], str)
            or not isinstance(term["relations"], list)
            or any(not isinstance(r, str) or not r.strip() for r in term["relations"])
        ):
            raise ValueError("Invalid saved GO term")
        ids.add(term["id"])
    if section["status"] == "no_matches" and section["terms"]:
        raise ValueError("Empty GO status has terms")
    if section["status"] != "unavailable" and section["error"] is not None:
        raise ValueError("Completed GO context has an error")
