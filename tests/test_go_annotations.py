"""GO terms, explicit relations, pagination, and offline index presentation."""

import copy
from dataclasses import asdict

import httpx
import pytest

from ligandscout.api_http import Audit
from ligandscout.go_annotations import enrich_go, validate_saved_go
from ligandscout.models import SearchOutcome, source_result_from_dict
from ligandscout.report import render_index
from ligandscout.storage import save_results
from ligandscout.uniprot import fetch_uniprot, parse_entry


def entry_with_go(entry):
    entry = copy.deepcopy(entry)
    entry["uniProtKBCrossReferences"] = [
        {
            "database": "GO",
            "id": go_id,
            "properties": [{"key": "GoTerm", "value": term}],
        }
        for go_id, term in (
            ("GO:0005524", "F:ATP binding"),
            ("GO:0006979", "P:response to oxidative stress"),
            ("GO:0005634", "C:nucleus"),
        )
    ]
    return entry


def annotation(go_id="GO:0005524", qualifier="enables", aspect="molecular_function"):
    return {
        "geneProductId": "UniProtKB:P00519",
        "goId": go_id,
        "qualifier": qualifier,
        "goAspect": aspect,
    }


def enrich(positive, uniprot_entry, now, handler):
    parsed = parse_entry(entry_with_go(uniprot_entry), query="P00519")
    # Use a normal adapter result to retain schema metadata.
    with httpx.Client(
        transport=httpx.MockTransport(lambda _: httpx.Response(200, json=uniprot_entry))
    ) as client:
        out = fetch_uniprot(
            "P00519", client=client, sleep=lambda _: None, now=now, monotonic=lambda: 0
        )
    out.result.source_context = parsed.source_context
    audit = Audit()
    with httpx.Client(transport=httpx.MockTransport(handler)) as client:
        enrich_go(
            out.result,
            entry_with_go(uniprot_entry),
            audit=audit,
            client=client,
            sleep=lambda _: None,
            now=now,
            monotonic=lambda: 0,
        )
    return SearchOutcome(out.result, None, audit.exchanges)


def test_paginated_relations_and_one_line_per_term(
    positive, uniprot_entry, now, tmp_path
):
    pages = []

    def handler(request):
        page = int(request.url.params["page"])
        pages.append(page)
        assert request.url.params["geneProductId"] == "UniProtKB:P00519"
        rows = (
            [annotation(), annotation(qualifier="contributes_to")]
            if page == 1
            else [annotation(qualifier="NOT|enables")]
        )
        return httpx.Response(
            200,
            json={
                "numberOfHits": 3,
                "results": rows,
                "pageInfo": {"current": page, "total": 2},
            },
        )

    out = enrich(positive, uniprot_entry, now, handler)
    assert pages == [1, 2]
    section = out.result.source_context["go_annotations"]
    validate_saved_go(section)
    assert section["status"] == "success"
    assert section["terms"][2]["relations"] == [
        "NOT|enables",
        "contributes_to",
        "enables",
    ]
    page = render_index(out)
    assert all(
        title in page
        for title in (
            "MF — Molecular Function",
            "BP — Biological Process",
            "CC — Cellular Component",
        )
    )
    assert page.count(">GO:0005524</a>") == 1
    assert "NOT|enables; contributes_to; enables" in page
    assert ">not supplied</td>" in page
    folder = save_results(out, root=tmp_path)
    assert (
        source_result_from_dict(asdict(out.result)).source_context
        == out.result.source_context
    )
    assert (folder / "sources/uniprot-requests.json").is_file()


@pytest.mark.parametrize(
    "problem", ["outage", "short", "identity", "aspect", "malformed"]
)
def test_optional_lookup_failure_keeps_terms(positive, uniprot_entry, now, problem):
    def handler(_):
        row = annotation()
        if problem == "outage":
            return httpx.Response(503)
        if problem == "identity":
            row["geneProductId"] = "UniProtKB:P00699"
        if problem == "aspect":
            row["goAspect"] = "biological_process"
        payload = {
            "numberOfHits": 2 if problem == "short" else 1,
            "results": [row],
            "pageInfo": {"current": 1, "total": 1},
        }
        return httpx.Response(200, json=[] if problem == "malformed" else payload)

    out = enrich(positive, uniprot_entry, now, handler)
    section = out.result.source_context["go_annotations"]
    assert section["status"] == "unavailable"
    assert len(section["terms"]) == 3
    assert all(not t["relations"] for t in section["terms"])
    assert "GO relation lookup unavailable" in render_index(out)
    assert out.result.error is None


def test_no_terms_avoids_quickgo(uniprot_entry, now):
    def handler(request):
        assert request.url.host == "rest.uniprot.org"
        data = copy.deepcopy(uniprot_entry)
        data.update(features=[], comments=[])
        return httpx.Response(200, json=data)

    with httpx.Client(transport=httpx.MockTransport(handler)) as client:
        out = fetch_uniprot(
            "P00519", client=client, sleep=lambda _: None, now=now, monotonic=lambda: 0
        )
    assert out.result.source_context["go_annotations"]["status"] == "no_matches"
    assert "No GO terms listed" in render_index(out)


def test_go_html_escaping_and_saved_validation(positive, uniprot_entry, now):
    out = enrich(
        positive,
        uniprot_entry,
        now,
        lambda _: httpx.Response(
            200,
            json={
                "numberOfHits": 0,
                "results": [],
                "pageInfo": {"current": 1, "total": 0},
            },
        ),
    )
    section = out.result.source_context["go_annotations"]
    section["terms"][0]["name"] = '<script>alert("x")</script>'
    assert "<script>" not in render_index(out)
    assert "&lt;script&gt;" in render_index(out)
    section["terms"][0]["id"] = '"><script>'
    with pytest.raises(ValueError):
        source_result_from_dict(asdict(out.result))


def test_fetch_pipeline_saves_go_and_reuses_without_network(
    uniprot_entry, now, tmp_path
):
    from ligandscout.cache import load_completed

    data = entry_with_go(uniprot_entry)
    data.update(features=[], comments=[])
    requests = []

    def handler(request):
        requests.append(request.url.host)
        if request.url.host == "rest.uniprot.org":
            return httpx.Response(200, json=data)
        assert request.url.host == "www.ebi.ac.uk"
        return httpx.Response(
            200,
            json={
                "numberOfHits": 1,
                "results": [annotation(qualifier="NOT|enables")],
                "pageInfo": {"current": 1, "total": 1},
            },
        )

    with httpx.Client(transport=httpx.MockTransport(handler)) as client:
        out = fetch_uniprot(
            "P00519", client=client, sleep=lambda _: None, now=now, monotonic=lambda: 0
        )
    assert requests == ["rest.uniprot.org", "www.ebi.ac.uk"]
    folder = save_results(out, root=tmp_path)
    reused = load_completed(tmp_path, "P00519", "UniProt")
    assert reused is not None and reused.reused
    assert "NOT|enables" in render_index(reused)
    assert ">not supplied</td>" in render_index(reused)
    assert "go_annotations" in (folder / "sources/uniprot-requests.json").read_text()
    # Legacy snapshots trigger one UniProt refresh to populate the new section.
    import json

    path = folder / "sources/uniprot.json"
    saved = json.loads(path.read_text())
    del saved["source_context"]["go_annotations"]
    path.write_text(json.dumps(saved))
    assert load_completed(tmp_path, "P00519", "UniProt") is None


def test_go_tables_have_two_aligned_columns(positive, uniprot_entry, now):
    from html.parser import HTMLParser

    class Tables(HTMLParser):
        def __init__(self):
            super().__init__()
            self.in_go = False
            self.rows = []

        def handle_starttag(self, tag, attrs):
            if tag == "table":
                self.in_go = dict(attrs).get("class") == "go-table"
            elif self.in_go and tag == "tr":
                self.rows.append([])
            elif self.in_go and tag in ("td", "th"):
                self.rows[-1].append(tag)

        def handle_endtag(self, tag):
            if tag == "table":
                self.in_go = False

    out = enrich(
        positive,
        uniprot_entry,
        now,
        lambda _: httpx.Response(
            200,
            json={
                "numberOfHits": 1,
                "results": [annotation()],
                "pageInfo": {"current": 1, "total": 1},
            },
        ),
    )
    page = render_index(out)
    parser = Tables()
    parser.feed(page)
    assert page.count('class="go-table"') == 3
    assert len(parser.rows) == 6  # header and one term for each aspect
    assert all(len(row) == 2 for row in parser.rows)
    assert "<td>enables</td>" in page
    assert page.count("width:72%") == 3


def test_index_function_description_source_and_escaping(uniprot_entry, now):
    data = copy.deepcopy(uniprot_entry)
    data.update(
        features=[],
        comments=[
            {
                "commentType": "FUNCTION",
                "texts": [
                    {"value": "Accession-level function <script>alert(1)</script>"}
                ],
            },
            {
                "commentType": "FUNCTION",
                "molecule": "Isoform <1>",
                "texts": [{"value": "Isoform-specific function."}],
            },
            {
                "commentType": "SUBUNIT",
                "texts": [{"value": "Subunit description should not be a function."}],
            },
        ],
    )
    with httpx.Client(
        transport=httpx.MockTransport(lambda _: httpx.Response(200, json=data))
    ) as client:
        out = fetch_uniprot(
            "P00519", client=client, sleep=lambda _: None, now=now, monotonic=lambda: 0
        )
    page = render_index(out)
    assert "Accession-level function &lt;script&gt;" in page
    assert "<script>" not in page
    assert "Isoform &lt;1&gt;" in page
    assert "Isoform-specific function." in page
    assert "Subunit description should not be a function." not in page
    assert "https://www.uniprot.org/uniprotkb/P00519/entry#function" in page
    assert "Source:" in page and "UniProtKB/Swiss-Prot (P00519)" in page
    assert page.index('id="protein-function"') < page.index('id="go-annotations"')


def test_index_function_missing_is_not_invented(uniprot_entry, now):
    data = copy.deepcopy(uniprot_entry)
    data.update(features=[], comments=[])
    with httpx.Client(
        transport=httpx.MockTransport(lambda _: httpx.Response(200, json=data))
    ) as client:
        out = fetch_uniprot(
            "P00519", client=client, sleep=lambda _: None, now=now, monotonic=lambda: 0
        )
    assert "No function description supplied" in render_index(out)
    out.result.error = type("Error", (), {"message": "UniProt unavailable"})()
    assert "Function description unavailable" in render_index(out)


def test_function_preserves_source_paragraphs_and_adds_spacing(uniprot_entry, now):
    from ligandscout.report import render_function

    data = copy.deepcopy(uniprot_entry)
    data.update(
        features=[],
        comments=[
            {
                "commentType": "FUNCTION",
                "texts": [
                    {
                        "value": "First source paragraph.\r\n\r\nSecond source paragraph.\nThird source paragraph."
                    },
                    {"value": "Separate source text <tag>."},
                ],
            }
        ],
    )
    with httpx.Client(
        transport=httpx.MockTransport(lambda _: httpx.Response(200, json=data))
    ) as client:
        out = fetch_uniprot(
            "P00519", client=client, sleep=lambda _: None, now=now, monotonic=lambda: 0
        )
    section = render_function(out.result)
    assert section.count('class="function-paragraph"') == 4
    assert (
        'First source paragraph.</p>\n\n<p class="function-paragraph">Second source paragraph.'
        in section
    )
    assert "Separate source text &lt;tag&gt;." in section
    page = render_index(out)
    assert ".function-paragraph{line-height:1.75;margin:0 0 1.75em}" in page
    assert "max-width:78ch" in page


def test_long_function_text_gets_reading_breaks_without_loss():
    from ligandscout.report import function_paragraphs

    text = (
        "A long function sentence "
        + "with source detail " * 40
        + ". Next complete sentence. Final source sentence."
    )
    paragraphs = list(function_paragraphs(text))
    assert len(paragraphs) == 2
    assert " ".join(paragraphs) == text
    assert all(paragraph.strip() for paragraph in paragraphs)
    # Source paragraph boundaries and abbreviation qualifiers stay intact.
    with_abbreviation = "Source detail " * 60 + "e.g. ATP binding. Another sentence."
    paragraphs = list(function_paragraphs(with_abbreviation))
    assert paragraphs[0].endswith("ATP binding.")
    assert " ".join(paragraphs) == with_abbreviation
    assert list(function_paragraphs("Short paragraph.\n\nNext paragraph.")) == [
        "Short paragraph.",
        "Next paragraph.",
    ]


@pytest.mark.parametrize(
    "entry_type,label",
    [
        ("UniProtKB reviewed (Swiss-Prot)", "UniProt | SwissProt"),
        (
            "UniProtKB unreviewed (TrEMBL)",
            'UniProt | <s title="Unreviewed UniProtKB/TrEMBL entry">SwissProt</s>',
        ),
        (None, "UniProt | SwissProt status unavailable"),
        ("Unexpected entry type", "UniProt | SwissProt status unavailable"),
    ],
)
def test_go_heading_shows_entry_review_status(entry_type, label):
    from types import SimpleNamespace
    from ligandscout.go_annotations import render_go

    result = SimpleNamespace(
        retrieved_at="2026-11-06T12:00:00Z",
        source_context={
            "entry": {"entry_type": entry_type},
            "go_annotations": {
                "status": "no_matches",
                "terms": [],
                "request_ids": [],
                "error": None,
            },
        },
    )
    page = render_go(result)
    assert f"<h2>GO annotations ({label})</h2>" in page
    assert "not the evidence level of each GO annotation" in page
    if entry_type != "UniProtKB unreviewed (TrEMBL)":
        assert "<s " not in page


def test_go_heading_absent_source_does_not_claim_swissprot():
    from ligandscout.go_annotations import render_go

    page = render_go(None)
    assert "SwissProt status unavailable" in page
    assert "<s " not in page
