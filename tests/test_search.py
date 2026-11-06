"""Ensure searches query all three implemented sources."""

import httpx

from ligandscout.search import run_search


def test_search_requests_all_sources(payload, now, empty_uniprot):
    requests = []

    def handler(request):
        requests.append(request)
        if request.url.host == "search.rcsb.org":
            return httpx.Response(204)
        if request.url.host == "rest.uniprot.org":
            return httpx.Response(200, json=empty_uniprot)
        assert request.url.host == "www.bindingdb.org"
        assert request.url.path == "/rest/getLigandsByUniprots"
        assert request.url.params["uniprot"] == "P00699"
        return httpx.Response(200, json=payload)

    with httpx.Client(transport=httpx.MockTransport(handler)) as client:
        outcome = run_search(" p00699 ", client=client, sleep=lambda _: None, now=now)
    assert len(requests) == 3
    assert outcome.accession == "P00699"
    assert [s.result.source for s in outcome.sources] == [
        "BindingDB",
        "UniProt",
        "RCSB PDB",
    ]
    assert outcome.sources[0].result.schema_version == 5
    assert outcome.sources[0].result.status == "success"
