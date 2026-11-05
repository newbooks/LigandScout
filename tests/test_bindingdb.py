"""Validate the API contract, evidence preservation, and failure handling."""

import copy

import httpx
import pytest

from ligandscout.accessions import normalize_accession
from ligandscout.bindingdb import URL


@pytest.mark.parametrize(
    "value,expected",
    [
        ("P00699", "P00699"),
        (" p00699 ", "P00699"),
        ("A0A024RBG1", "A0A024RBG1"),
        ("Q9Y261", "Q9Y261"),
    ],
)
def test_normalization(value, expected):
    assert normalize_accession(value) == expected


@pytest.mark.parametrize(
    "value",
    ["", "LYZ", "P00699-2", "../P00699", "P00699,P00176", "P006999", ">P00699\nABC"],
)
def test_invalid_accessions(value):
    with pytest.raises(ValueError):
        normalize_accession(value)


def test_request_and_provenance(fetch, payload):
    def handler(request):
        assert str(request.url).split("?")[0] == URL
        assert dict(request.url.params) == {
            "uniprot": "P00699",
            "response": "application/json",
        }
        assert request.headers["Accept"] == "application/json"
        assert request.headers["User-Agent"] == "ligandscout/0.1.0"
        assert request.extensions["timeout"] == dict.fromkeys(
            ["connect", "read", "write", "pool"], 30.0
        )
        return httpx.Response(200, json=payload, headers={"Content-Type": "text/html"})

    outcome, sleeps = fetch(handler)
    result = outcome.result
    assert result.status == "success"
    assert result.retrieved_at == "2026-11-04T15:00:00.000000Z"
    assert result.source_version is None and result.error is None
    record = result.records[0]
    assert record.compound_id == "50027784"
    assert record.measurement_value == "140" and record.measurement_unit is None
    assert record.pmid == "16248836" and record.doi.startswith("10.")
    assert record.compound_name is None and record.inchikey is None
    assert outcome.raw_body and not sleeps


@pytest.mark.parametrize(
    "body",
    [b'{"getLindsByUniprotsResponse":{"affinities":[]}}'],
)
def test_empty_variants(fetch, body):
    outcome, sleeps = fetch(lambda request: httpx.Response(200, content=body))
    assert outcome.result.status == "no_matches"
    assert outcome.result.records == [] and outcome.result.error is None
    assert outcome.raw_body == body and sleeps == []


@pytest.mark.parametrize(
    "body",
    [
        b"<html>unavailable</html>",
        b"{",
        b"null",
        b"[]",
        b"{}",
        b'{"getLindsByUniprotsResponse":null}',
        b'{"getLindsByUniprotsResponse":{"affinities":null}}',
        b'{"getLindsByUniprotsResponse":{"affinities":[null]}}',
        b'{"getLindsByUniprotsResponse":{"affinities":[]},"extra":NaN}',
        b'"   "',
    ],
)
def test_schema_failure(fetch, body):
    outcome, sleeps = fetch(lambda request: httpx.Response(200, content=body))
    assert outcome.result.status == "failed"
    assert outcome.result.error.kind == "schema"
    assert outcome.result.error.http_status == 200
    assert outcome.result.records == [] and outcome.raw_body == body and not sleeps
    assert "<html>" not in outcome.result.error.message


@pytest.mark.parametrize(
    "field,value",
    [
        ("monomerid", 1),
        ("monomerid", "１２"),
        ("monomerid", "../1"),
        ("affinity_type", " "),
        ("affinity", None),
        ("affinity", 10),
        ("pmid", []),
        ("doi", 1),
        ("query", {}),
        ("smile", False),
    ],
)
def test_invalid_record_fails_entire_source(fetch, payload, field, value):
    row = copy.deepcopy(payload["getLindsByUniprotsResponse"]["affinities"][0])
    row[field] = value
    payload["getLindsByUniprotsResponse"]["affinities"].append(row)
    outcome, _ = fetch(lambda request: httpx.Response(200, json=payload))
    assert outcome.result.status == "failed" and outcome.result.records == []


def test_preserves_qualifiers_duplicates_unknown_types_and_optional_fields(
    fetch, payload
):
    row = {
        "monomerid": "2",
        "affinity_type": "Novel",
        "affinity": "> 1",
        "pmid": " ",
        "extra": 1,
    }
    payload["getLindsByUniprotsResponse"]["affinities"] += [row, row]
    outcome, _ = fetch(lambda request: httpx.Response(200, json=payload))
    assert len(outcome.result.records) == 3
    assert outcome.result.records[0].measurement_value == "> 1"
    assert outcome.result.records[0].smiles is None
    assert outcome.result.records[0].pmid is None
    assert outcome.result.records[0] == outcome.result.records[1]


@pytest.mark.parametrize(
    "failure",
    [httpx.ConnectError, httpx.ReadError, httpx.ReadTimeout, httpx.PoolTimeout],
)
def test_transport_exhaustion(fetch, failure):
    calls = []

    def handler(request):
        calls.append(request)
        raise failure("secret payload must not be logged", request=request)

    outcome, sleeps = fetch(handler)
    assert len(calls) == 3 and sleeps == [1, 2]
    assert outcome.result.error.kind == "transport"
    assert outcome.result.error.http_status is None and outcome.raw_body is None
    assert "secret" not in outcome.result.error.message


def test_protocol_error_not_retried(fetch):
    outcome, sleeps = fetch(
        lambda request: (_ for _ in ()).throw(httpx.RemoteProtocolError("bad"))
    )
    assert outcome.result.error.kind == "transport" and sleeps == []


@pytest.mark.parametrize("status", [400, 404, 204])
def test_nonretryable_http_failure_even_when_empty(fetch, status):
    outcome, sleeps = fetch(lambda request: httpx.Response(status, content=b""))
    assert outcome.result.error.kind == "http"
    assert outcome.result.error.http_status == status and sleeps == []


@pytest.mark.parametrize("status", [429, 500, 503, 599])
def test_http_retry_exhaustion(fetch, status):
    outcome, sleeps = fetch(lambda request: httpx.Response(status, content=b"down"))
    assert sleeps == [1, 2]
    assert outcome.result.error.http_status == status and outcome.raw_body == b"down"


@pytest.mark.parametrize(
    "header,expected",
    [
        ("4", [4, 4]),
        ("0", [1, 2]),
        ("nonsense", [1, 2]),
        ("-1", [1, 2]),
        ("Wed, 04 Nov 2026 15:00:03 GMT", [3, 3]),
        ("Wed, 04 Nov 2026 14:59:00 GMT", [1, 2]),
        ("61", []),
    ],
)
def test_retry_after(fetch, header, expected):
    outcome, sleeps = fetch(
        lambda request: httpx.Response(429, headers={"Retry-After": header})
    )
    assert sleeps == expected
    if header == "61":
        assert "Retry later" in outcome.result.error.message


def test_retry_success_clears_previous_error(fetch, payload):
    calls = []

    def handler(request):
        calls.append(request)
        return (
            httpx.Response(503)
            if len(calls) == 1
            else httpx.Response(200, json=payload)
        )

    outcome, sleeps = fetch(handler)
    assert outcome.result.status == "success" and outcome.result.error is None
    assert sleeps == [1]


@pytest.mark.parametrize(
    "location",
    [
        "http://example.org/unsafe",
        "https://",
        "https://[",
        "https://user:pass@example.org/",
        None,
    ],
)
def test_invalid_redirect(fetch, location):
    headers = {"Location": location} if location is not None else {}
    outcome, sleeps = fetch(
        lambda request: httpx.Response(302, headers=headers, content=b"redirect")
    )
    assert outcome.result.error.kind == "http" and not sleeps
    assert outcome.raw_body == b"redirect"


def test_relative_https_redirect(fetch, payload):
    def handler(request):
        if request.url.path == "/rest/getLigandsByUniprots":
            return httpx.Response(302, headers={"Location": "/new"})
        assert request.url.scheme == "https" and request.url.path == "/new"
        return httpx.Response(200, json=payload)

    outcome, sleeps = fetch(handler)
    assert outcome.result.status == "success" and not sleeps
    assert outcome.result.request["url"] == URL


def test_redirect_loop_stops_and_does_not_retry(fetch):
    calls = []

    def handler(request):
        calls.append(request)
        return httpx.Response(302, headers={"Location": "/loop"})

    outcome, sleeps = fetch(handler)
    assert len(calls) == 6 and not sleeps
    assert outcome.result.error.kind == "http"


def test_redirect_does_not_reset_retries(fetch):
    paths = []

    def handler(request):
        paths.append(request.url.path)
        if request.url.path == "/rest/getLigandsByUniprots":
            return httpx.Response(302, headers={"Location": "/new"})
        return httpx.Response(503)

    outcome, sleeps = fetch(handler)
    assert len(paths) == 6 and sleeps == [1, 2]
    assert outcome.result.error.http_status == 503
    assert paths[::2] == ["/rest/getLigandsByUniprots"] * 3


def test_transport_failure_discards_earlier_body(fetch):
    calls = []

    def handler(request):
        calls.append(request)
        if len(calls) == 1:
            return httpx.Response(503, content=b"old")
        raise httpx.ReadTimeout("timeout")

    outcome, sleeps = fetch(handler)
    assert outcome.raw_body is None and outcome.result.error.kind == "transport"
    assert sleeps == [1, 2]


@pytest.mark.parametrize("body", [b"", b" \r\n\t", b'""'])
def test_ambiguous_empty_retries_then_fails(fetch, body):
    outcome, sleeps = fetch(lambda request: httpx.Response(200, content=body))
    assert sleeps == [1, 2]
    assert outcome.result.status == "failed"
    assert outcome.result.error.kind == "schema"
    assert "cannot determine whether matches exist" in outcome.result.error.message
    assert outcome.result.records == [] and outcome.raw_body == body


@pytest.mark.parametrize("final", ["success", "no_matches"])
def test_empty_then_structured_result_recovers(fetch, payload, final):
    if final == "no_matches":
        payload["getLindsByUniprotsResponse"]["affinities"] = []
    calls = []

    def handler(request):
        calls.append(request)
        return (
            httpx.Response(200, content=b"")
            if len(calls) == 1
            else httpx.Response(200, json=payload)
        )

    outcome, sleeps = fetch(handler)
    assert len(calls) == 2 and sleeps == [1]
    assert outcome.result.status == final and outcome.result.error is None


def test_http_and_empty_responses_share_attempt_budget(fetch):
    calls = []

    def handler(request):
        calls.append(request)
        return (
            httpx.Response(503) if len(calls) == 1 else httpx.Response(200, content=b"")
        )

    outcome, sleeps = fetch(handler)
    assert len(calls) == 3 and sleeps == [1, 2]
    assert outcome.result.error.kind == "schema"
