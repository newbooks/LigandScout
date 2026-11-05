"""Ensure simple mode queries only the implemented source."""

import httpx

from ligandscout.search import run_simple


def test_simple_mode_requests_bindingdb_and_uniprot(payload, now, empty_uniprot):
    requests = []

    def handler(request):
        requests.append(request)
        if request.url.host == "rest.uniprot.org":
            return httpx.Response(200, json=empty_uniprot)
        assert request.url.host == "www.bindingdb.org"
        assert request.url.path == "/rest/getLigandsByUniprots"
        assert request.url.params["uniprot"] == "P00699"
        return httpx.Response(200, json=payload)

    with httpx.Client(transport=httpx.MockTransport(handler)) as client:
        outcome = run_simple(" p00699 ", client=client, sleep=lambda _: None, now=now)
    assert len(requests) == 2
    assert outcome.accession == "P00699"
    assert [s.result.source for s in outcome.sources] == ["BindingDB", "UniProt"]
    assert outcome.sources[0].result.schema_version == 4
    assert outcome.sources[0].result.status == "success"
