"""Offline ChEMBL API, selection, provenance and failure regression tests."""

import copy
import csv
import hashlib
import io
import json
from dataclasses import asdict, replace
from pathlib import Path

import httpx
import pytest
from typer.testing import CliRunner

from ligandscout.chembl import fetch_chembl
from ligandscout.models import RunOutcome, source_result_from_dict, upgrade_bindingdb
from ligandscout.report import render_html, render_summary
from ligandscout.storage import save_results


class Clock:
    def __init__(self):
        self.time = 0.0
        self.sleeps = []

    def monotonic(self):
        return self.time

    def sleep(self, delay):
        self.sleeps.append(delay)
        self.time += delay


class API:
    def __init__(self, data, clock, page_size=2):
        self.data, self.clock, self.page_size = data, clock, page_size
        self.requests, self.starts = [], []
        self.versions = ["ChEMBL_test", "ChEMBL_test"]
        self.status_calls = 0
        self.accession = "P00533"
        self.hook = None

    def __call__(self, request):
        self.requests.append(request)
        self.starts.append(self.clock.time)
        if self.hook:
            response = self.hook(request)
            if response is not None:
                return response
        return self.response(request)

    def response(self, request):
        resource = request.url.path.rsplit("/", 1)[-1].removesuffix(".json")
        params = request.url.params
        if resource == "status":
            v = self.versions[min(self.status_calls, 1)]
            self.status_calls += 1
            return httpx.Response(200, json={"status": "UP", "chembl_db_version": v})
        key = "activities" if resource == "activity" else resource + "s"
        rows = self.data[key]
        id_field = {"activity": "activity_id"}.get(resource, resource + "_chembl_id")
        assert params["order_by"] == id_field
        if resource == "target":
            assert params["target_components__accession"] == self.accession
        elif resource == "activity":
            rows = [
                r for r in rows if r["target_chembl_id"] == params["target_chembl_id"]
            ]
        else:
            ids = params[id_field + "__in"].split(",")
            assert len(ids) <= 50
            rows = [r for r in rows if r[id_field] in ids]
        offset = int(params["offset"])
        limit = min(int(params["limit"]), self.page_size)
        url = request.url.copy_set_param("offset", str(offset + limit)).copy_set_param(
            "limit", str(limit)
        )
        return httpx.Response(
            200,
            json={
                key: rows[offset : offset + limit],
                "page_meta": {
                    "offset": offset,
                    "limit": limit,
                    "total_count": len(rows),
                    "next": str(url) if offset + limit < len(rows) else None,
                },
            },
        )


@pytest.fixture
def chembl_data():
    return json.loads(
        (Path(__file__).parent / "fixtures/chembl-search.json").read_text()
    )


@pytest.fixture
def chembl_case(chembl_data, now):
    def call(configure=None, accession="P00533", data=None):
        clock = Clock()
        api = API(copy.deepcopy(chembl_data if data is None else data), clock)
        api.accession = accession
        if configure:
            configure(api)
        with httpx.Client(transport=httpx.MockTransport(api)) as client:
            outcome = fetch_chembl(
                accession,
                client=client,
                sleep=clock.sleep,
                now=now,
                monotonic=clock.monotonic,
            )
            assert not client.is_closed
        return outcome, api, clock

    return call


def test_acceptance_selection_and_pagination(chembl_case):
    outcome, api, clock = chembl_case()
    r = outcome.result
    assert r.status == "success", r.error
    assert r.counts == {
        "complete": True,
        "candidate_targets": 2,
        "eligible_targets": 1,
        "candidate_activities": 5,
        "included_activities": 3,
        "excluded_activities": 2,
        "unverified_activities": 1,
        "unique_ligands": 2,
    }
    assert [row.source_record_id for row in r.records] == ["1", "2", "5"]
    audit = {row["activity_id"]: row for row in r.selection_audit["activities"]}
    assert audit["3"]["reasons"] == ["confidence_not_9", "homolog_assignment"]
    assert audit["4"]["reasons"] == ["unverified_target_assignment"]
    assert r.selection_audit["targets"][1]["reasons"] == [
        "component_count",
        "target_type",
    ]
    assert not any(
        "CHEMBL2111431" in str(q.url) and q.url.path.endswith("activity.json")
        for q in api.requests
    )
    assert min(b - a for a, b in zip(api.starts, api.starts[1:])) >= 0.5
    assert clock.sleeps
    assert r.source_version == "ChEMBL_test"
    assert len(outcome.exchanges) == len(api.requests)
    assert r.records[0].measurement_value == "0"
    assert r.records[0].measurement_relation == ">"
    assert r.records[0].source_values["published"]["units"] == "uM"
    assert (
        r.records[1].variant_mutation == "L858R" and r.records[1].variant_isoform == "1"
    )
    assert r.records[2].measurement_basis == "qualitative"
    assert r.records[2].measurement_text == "Inactive"
    assert r.records[2].measurement_value is None
    assert r.records[0].smiles is None and r.records[0].compound_name is None
    assert r.records[0].pmid == "12345"


@pytest.mark.parametrize(
    "mode,reason",
    [
        ("targets", "no_targets"),
        ("eligible", "no_eligible_targets"),
        ("activities", "no_activities"),
        ("excluded", "all_activities_excluded"),
    ],
)
def test_no_matches_cases(chembl_case, mode, reason):
    def configure(api):
        if mode == "targets":
            api.data["targets"] = []
        elif mode == "eligible":
            api.data["targets"][0]["target_type"] = "PROTEIN FAMILY"
        elif mode == "activities":
            api.data["activities"] = []
        else:
            for a in api.data["assays"]:
                a["relationship_type"] = "H"

    outcome, _, _ = chembl_case(configure)
    assert outcome.result.status == "no_matches"
    assert outcome.result.selection_audit["outcome_reason"] == reason
    assert outcome.result.counts["unique_ligands"] == 0


def test_no_related_accession_fallback(chembl_case):
    def configure(api):
        api.hook = lambda r: (
            httpx.Response(
                200,
                json={
                    "targets": [],
                    "page_meta": {
                        "offset": 0,
                        "limit": 100,
                        "total_count": 0,
                        "next": None,
                    },
                },
            )
            if r.url.path.endswith("target.json")
            else None
        )

    outcome, api, _ = chembl_case(configure, accession="P0DTC1")
    assert outcome.result.status == "no_matches"
    assert all("P0DTD1" not in str(r.url) for r in api.requests)


@pytest.mark.parametrize(
    "change",
    [
        lambda d: d["activities"][0].update(activity_id=True),
        lambda d: d["activities"][0].update(standard_value=[]),
        lambda d: d["activities"][2].update(
            potential_duplicate=2
        ),  # excluded row still validated
        lambda d: d["assays"][0].update(confidence_score=True),
        lambda d: d["assays"][0].update(confidence_score=10),
        lambda d: d["targets"][0].update(target_components=None),
        lambda d: d["assays"][0].update(target_chembl_id="CHEMBL999"),
        lambda d: d["assays"][0].update(assay_type="F"),
        lambda d: d["activities"][1].update(assay_variant_mutation="T790M"),
        lambda d: d["molecules"].pop(),
        lambda d: d["documents"].clear(),
        lambda d: d["assays"].pop(),
        lambda d: d["molecules"][0].update(molecule_structures=[]),
    ],
)
def test_bad_data_is_failed_not_empty_or_partial(chembl_case, change):
    def configure(api):
        change(api.data)

    outcome, api, _ = chembl_case(configure)
    assert outcome.result.status == "failed"
    assert outcome.result.error.kind == "schema"
    assert outcome.result.records == []
    assert outcome.result.counts["unique_ligands"] is None
    assert outcome.result.selection_audit["outcome_reason"] == "incomplete"
    assert api.status_calls == 2


@pytest.mark.parametrize(
    "field,value",
    [
        ("relationship_type", None),
        ("confidence_score", None),
        ("assay_type", None),
        ("relationship_type", "UNKNOWN"),
        ("confidence_score", 0),
        ("assay_type", "A"),
    ],
)
def test_unverified_and_unsupported_assignments_excluded(chembl_case, field, value):
    def configure(api):
        api.data["assays"][0][field] = value
        if field == "assay_type":
            api.data["activities"][0]["assay_type"] = value

    outcome, _, _ = chembl_case(configure)
    assert outcome.result.status == "success"
    assert "1" not in {r.source_record_id for r in outcome.result.records}


@pytest.mark.parametrize(
    "bad_next",
    [
        "http://www.ebi.ac.uk/chembl/api/data/activity.json?offset=2&limit=2",
        "https://evil.invalid/chembl/api/data/activity.json?offset=2&limit=2",
        "https://www.ebi.ac.uk/chembl/api/data/molecule.json?offset=2&limit=2",
        "?offset=0&limit=2&target_chembl_id=CHEMBL203&order_by=activity_id",
        "?offset=2&limit=2&target_chembl_id=CHEMBL999&order_by=activity_id",
        "?offset=2&limit=2&target_chembl_id=CHEMBL203&order_by=activity_id&limit=2",
        "?offset=2&limit=2&target_chembl_id=CHEMBL203&order_by=activity_id#fragment",
        "?offset=99&limit=2&target_chembl_id=CHEMBL203&order_by=activity_id",
    ],
)
def test_unsafe_pagination_fails(chembl_case, bad_next):
    def configure(api):
        def hook(request):
            if request.url.path.endswith("activity.json"):
                d = json.loads(api.response(request).content)
                d["page_meta"]["next"] = bad_next
                return httpx.Response(200, json=d)

        api.hook = hook

    outcome, _, _ = chembl_case(configure)
    assert outcome.result.status == "failed"
    assert outcome.result.error.kind == "schema"


@pytest.mark.parametrize(
    "mode", ["total", "empty", "truncated", "conflict", "offset", "boolean"]
)
def test_incomplete_pagination(chembl_case, mode):
    def configure(api):
        def hook(request):
            if request.url.path.endswith("activity.json"):
                d = json.loads(api.response(request).content)
                offset = int(request.url.params["offset"])
                if mode == "total" and offset:
                    d["page_meta"]["total_count"] += 1
                elif mode == "empty":
                    d["activities"] = []
                elif mode == "truncated":
                    d["page_meta"]["next"] = None
                elif mode == "conflict" and offset:
                    d["activities"][0]["activity_id"] = 1
                    d["activities"][0]["standard_value"] = "different"
                elif mode == "offset":
                    d["page_meta"]["offset"] = 999
                elif mode == "boolean":
                    d["page_meta"]["limit"] = True
                return httpx.Response(200, json=d)

        api.hook = hook

    outcome, _, _ = chembl_case(configure)
    assert outcome.result.status == "failed"


@pytest.mark.parametrize(
    "body", [b"", b'""', b"<html>error</html>", b'{"value":NaN}', b'{"targets":{}}']
)
def test_bad_success_bodies_fail_without_retry(chembl_case, body):
    def configure(api):
        api.hook = lambda r: (
            httpx.Response(200, content=body)
            if r.url.path.endswith("target.json")
            else None
        )

    outcome, api, _ = chembl_case(configure)
    assert outcome.result.error.kind == "schema"
    assert len([r for r in api.requests if r.url.path.endswith("target.json")]) == 1


@pytest.mark.parametrize("failure", [429, 500, 503, "timeout", "network"])
def test_retries_and_preserves_all_exchanges(chembl_case, failure):
    def configure(api):
        attempts = []

        def hook(request):
            if request.url.path.endswith("target.json") and len(attempts) < 2:
                attempts.append(1)
                if failure == "timeout":
                    raise httpx.ReadTimeout("slow")
                if failure == "network":
                    raise httpx.ConnectError("offline")
                return httpx.Response(failure, content=b"error")

        api.hook = hook

    outcome, _, clock = chembl_case(configure)
    assert outcome.result.status == "success", outcome.result.error
    targets = [e for e in outcome.exchanges if e.role == "targets"]
    assert [e.attempt for e in targets] == [1, 2, 3]
    assert len({e.request_id for e in targets}) == 1
    assert 1 in clock.sleeps and 2 in clock.sleeps
    if isinstance(failure, str):
        assert targets[0].body is None and targets[0].transport_error
    else:
        assert targets[0].body == b"error"


@pytest.mark.parametrize("status", [400, 404, 429, 500])
def test_exhausted_or_nonretryable_http(chembl_case, status):
    def configure(api):
        api.hook = lambda r: (
            httpx.Response(status, content=b"error")
            if r.url.path.endswith("target.json")
            else None
        )

    outcome, api, _ = chembl_case(configure)
    assert outcome.result.error.kind == "http"
    assert outcome.result.error.http_status == status
    assert len([r for r in api.requests if r.url.path.endswith("target.json")]) == (
        3 if status in (429, 500) else 1
    )


@pytest.mark.parametrize("retry_after", ["61", "999999999999999999999"])
def test_retry_after_above_limit_stops(chembl_case, retry_after):
    def configure(api):
        api.hook = lambda r: (
            httpx.Response(429, headers={"Retry-After": retry_after})
            if r.url.path.endswith("target.json")
            else None
        )

    outcome, api, clock = chembl_case(configure)
    assert outcome.result.error.kind == "http"
    assert "60 seconds" in outcome.result.error.message
    assert all(delay <= 60 for delay in clock.sleeps)
    assert len([r for r in api.requests if r.url.path.endswith("target.json")]) == 1


@pytest.mark.parametrize(
    "header,expected",
    [("4", 4), ("Wed, 04 Nov 2026 15:00:07 GMT", 7), ("bad", 1), ("-3", 1)],
)
def test_retry_after_delay(chembl_case, header, expected):
    def configure(api):
        calls = []

        def hook(r):
            if r.url.path.endswith("target.json") and not calls:
                calls.append(1)
                return httpx.Response(429, headers={"Retry-After": header})

        api.hook = hook

    outcome, _, clock = chembl_case(configure)
    assert outcome.result.status == "success"
    assert expected in clock.sleeps


@pytest.mark.parametrize("change", ["version", "missing", "outage", "wrong_type"])
def test_release_metadata(chembl_case, change):
    def configure(api):
        if change == "version":
            api.versions[1] = "new_version"
        elif change == "missing":
            api.versions[1] = None
        elif change == "wrong_type":
            api.versions[0] = 37
        else:
            api.hook = lambda r: (
                httpx.Response(404) if r.url.path.endswith("status.json") else None
            )

    outcome, _, _ = chembl_case(configure)
    if change == "version":
        assert outcome.result.status == "failed"
        assert "release changed" in outcome.result.error.message
    else:
        assert outcome.result.status == "success"
        assert outcome.result.warnings
    assert outcome.result.source_version is None


@pytest.mark.parametrize(
    "location",
    [
        "http://www.ebi.ac.uk/chembl/api/data/target.json",
        "https://evil.invalid/target.json",
        None,
        "loop",
    ],
)
def test_invalid_or_looping_redirect(chembl_case, location):
    def configure(api):
        def hook(r):
            if r.url.path.endswith("target.json"):
                return httpx.Response(
                    302,
                    headers={}
                    if location is None
                    else {"Location": str(r.url) if location == "loop" else location},
                )

        api.hook = hook

    outcome, api, _ = chembl_case(configure)
    assert outcome.result.error.kind == "http"
    assert len([r for r in api.requests if r.url.path.endswith("target.json")]) == (
        6 if location == "loop" else 1
    )


def test_safe_redirect_and_retry_counter(chembl_case):
    def configure(api):
        calls = []

        def hook(r):
            if r.url.path.endswith("target.json"):
                calls.append(1)
                if len(calls) in (1, 3):
                    return httpx.Response(302, headers={"Location": str(r.url)})
                if len(calls) == 2:
                    return httpx.Response(503)

        api.hook = hook

    outcome, _, _ = chembl_case(configure)
    assert outcome.result.status == "success"
    assert [
        (e.attempt, e.redirect_hop) for e in outcome.exchanges if e.role == "targets"
    ] == [(1, 0), (1, 1), (2, 0), (2, 1)]


def test_batch_boundary_and_noncollapsed_parent(chembl_case, chembl_data):
    data = copy.deepcopy(chembl_data)
    activity = copy.deepcopy(data["activities"][0])
    data["activities"] = []
    data["molecules"] = []
    for i in range(51):
        data["activities"].append(
            dict(activity, activity_id=i + 1, molecule_chembl_id=f"CHEMBL{1000 + i}")
        )
        data["molecules"].append(
            {
                "molecule_chembl_id": f"CHEMBL{1000 + i}",
                "molecule_hierarchy": {"parent_chembl_id": "CHEMBL99"},
            }
        )
    outcome, api, _ = chembl_case(data=data)
    assert outcome.result.counts["unique_ligands"] == 51
    initial = [
        r
        for r in api.requests
        if r.url.path.endswith("molecule.json") and r.url.params["offset"] == "0"
    ]
    assert [
        len(r.url.params["molecule_chembl_id__in"].split(",")) for r in initial
    ] == [50, 1]


def test_publication_failure_preserves_evidence_and_other_source(
    chembl_case, positive, tmp_path
):
    outcome, _, _ = chembl_case()
    binding = replace(
        positive, result=replace(upgrade_bindingdb(positive.result), accession="P00533")
    )
    run = RunOutcome("P00533", [binding, outcome])
    folder = save_results(run, root=tmp_path)
    report = (folder / "report.md").read_text()
    assert "| ChEMBL | 2 | [View ligands](chembl-ligands.html) |" in report
    assert "| BindingDB | 1 |" in report
    page = (folder / "chembl-ligands.html").read_text()
    assert page.count("<tr><td><a ") == 2
    assert "/explore/compound/CHEMBL1" in page and "monomerid=" not in page
    assert "L858R" in page and "Inactive" in page and "&lt;Example ligand&gt;" in page
    manifest = json.loads((folder / "sources/chembl-requests.json").read_text())
    for entry, exchange in zip(manifest["exchanges"], outcome.exchanges):
        assert entry["url"] == exchange.url
        if exchange.body is not None:
            raw = (folder / "sources" / entry["raw_path"]).read_bytes()
            assert raw == exchange.body
            assert hashlib.sha256(raw).hexdigest() == entry["sha256"]
    rows = list(
        csv.DictReader(
            io.StringIO((folder / "ligands.tsv").read_text()), delimiter="\t"
        )
    )
    assert [r["source"] for r in rows] == ["BindingDB", "ChEMBL", "ChEMBL", "ChEMBL"]
    assert json.loads(rows[1]["target_accessions"]) == ["P00533"]
    assert rows[1]["potential_duplicate"] == "false"
    saved = json.loads((folder / "sources/chembl.json").read_text())
    assert asdict(source_result_from_dict(saved)) == saved
    assert "Excluded activities: 2" in render_summary(run, folder=tmp_path)


def test_failed_rerun_replaces_stale_html_and_raw(chembl_case, tmp_path):
    outcome, _, _ = chembl_case()
    folder = save_results(outcome, root=tmp_path)
    (folder / "sources/chembl-raw/keep.txt").write_text("user data")
    previous = set((folder / "sources/chembl-raw").glob("*.json"))

    def configure(api):
        api.hook = lambda r: (
            httpx.Response(404) if r.url.path.endswith("target.json") else None
        )

    failed, _, _ = chembl_case(configure)
    save_results(failed, root=tmp_path)
    assert "| ChEMBL | Unavailable |" in (folder / "report.md").read_text()
    assert "<tr><td><a " not in (folder / "chembl-ligands.html").read_text()
    assert (
        len(
            list(
                csv.DictReader(
                    io.StringIO((folder / "ligands.tsv").read_text()), delimiter="\t"
                )
            )
        )
        == 0
    )
    assert (folder / "sources/chembl-raw/keep.txt").read_text() == "user data"
    manifest = json.loads((folder / "sources/chembl-requests.json").read_text())
    expected = {
        folder / "sources" / e["raw_path"]
        for e in manifest["exchanges"]
        if e["raw_path"]
    }
    assert {
        p for p in (folder / "sources/chembl-raw").iterdir() if p.name != "keep.txt"
    } == expected
    assert previous - expected


def test_saved_v1_reader_and_invalid_versions(positive):
    data = asdict(positive.result)
    for key in ("policy", "counts", "warnings", "selection_audit"):
        data.pop(key)
    original = list(data["records"][0])[:12]
    data["records"] = [{k: r[k] for k in original} for r in data["records"]]
    r = source_result_from_dict(data)
    assert r.schema_version == 2 and r.counts["unique_ligands"] == 1
    assert r.records[0].target_match_status == "unverified"
    for key, value in [
        ("schema_version", 3),
        ("source", "ChEMBL"),
        ("schema_version", True),
    ]:
        with pytest.raises(ValueError):
            source_result_from_dict(dict(data, **{key: value}))


def test_multiple_eligible_targets_and_cross_target_join(chembl_case, chembl_data):
    data = copy.deepcopy(chembl_data)
    target = copy.deepcopy(data["targets"][0])
    target["target_chembl_id"] = "CHEMBL204"
    data["targets"].append(target)
    assay = dict(
        data["assays"][0], assay_chembl_id="CHEMBL506", target_chembl_id="CHEMBL204"
    )
    data["assays"].append(assay)
    data["activities"].append(
        dict(
            data["activities"][0],
            activity_id=6,
            assay_chembl_id="CHEMBL506",
            target_chembl_id="CHEMBL204",
        )
    )
    outcome, _api, _ = chembl_case(data=data)
    assert outcome.result.status == "success"
    assert outcome.result.counts["eligible_targets"] == 2
    assert outcome.result.counts["included_activities"] == 4
    assert outcome.result.counts["unique_ligands"] == 2
    assert {r.target_id for r in outcome.result.records} == {"CHEMBL203", "CHEMBL204"}


@pytest.mark.parametrize(
    "field,value",
    [
        ("component_type", "RNA"),
        ("accession", "P99999"),
        ("accession", None),
        ("component_type", None),
    ],
)
def test_component_identity_checks(chembl_case, field, value):
    def configure(api):
        api.data["targets"][0]["target_components"][0][field] = value

    outcome, _, _ = chembl_case(configure)
    assert outcome.result.status == "no_matches"
    assert outcome.result.selection_audit["outcome_reason"] == "no_eligible_targets"


def test_qualitative_fallback_unknown_endpoint_and_range(chembl_case):
    def configure(api):
        api.data["activities"][0].update(
            standard_value=None,
            type="Unknown endpoint",
            value=0,
            relation="<=",
            upper_value="2",
            units=None,
            text_value="range note",
        )
        api.data["activities"][1].update(standard_units=None, standard_upper_value="30")

    outcome, _, _ = chembl_case(configure)
    first, second, _ = outcome.result.records
    assert first.measurement_basis == "published" and first.measurement_value == "0"
    assert (
        first.measurement_type == "Unknown endpoint" and first.measurement_unit is None
    )
    assert (
        first.measurement_upper_value == "2" and first.measurement_text == "range note"
    )
    assert second.measurement_basis == "standard" and second.measurement_unit is None
    page = render_html(outcome.result)
    assert "&lt;= 0 to 2; range note" in page
    assert "12.5 to 30" in page


def test_lossless_json_numbers(chembl_case):
    def configure(api):
        def hook(request):
            if request.url.path.endswith("activity.json"):
                d = api.response(request)
                return httpx.Response(
                    200,
                    content=d.content.replace(
                        b'"standard_value":"0"',
                        b'"standard_value":0.12345678901234567890123456789',
                    ),
                )

        api.hook = hook

    outcome, _, _ = chembl_case(configure)
    assert outcome.result.status == "success"
    assert (
        outcome.result.records[0].measurement_value == "0.12345678901234567890123456789"
    )


def test_variant_activity_fallback_and_no_type_or_value(chembl_case):
    def configure(api):
        api.data["activities"][0]["assay_variant_mutation"] = "Deletion"
        api.data["activities"][4].update(
            type=None, standard_type=None, activity_comment=None
        )
        api.data["activities"][1].pop("assay_type")

    outcome, _, _ = chembl_case(configure)
    assert outcome.result.status == "success"
    assert outcome.result.records[0].variant_mutation == "Deletion"
    assert outcome.result.records[2].measurement_type is None
    assert "Unspecified" in render_summary(outcome, folder=Path("/tmp"))
    assert "None:" not in render_html(outcome.result)


def test_missing_optional_document_and_names(chembl_case):
    def configure(api):
        for a in api.data["activities"] + api.data["assays"]:
            a.pop("document_chembl_id", None)
        api.data["molecules"][0]["molecule_hierarchy"] = None

    outcome, api, _ = chembl_case(configure)
    assert outcome.result.status == "success"
    assert outcome.result.records[0].pmid is None
    assert not any(r.url.path.endswith("document.json") for r in api.requests)


@pytest.mark.parametrize("resource", ["assay", "molecule", "document"])
def test_unexpected_batch_entity_fails(chembl_case, resource):
    def configure(api):
        def hook(r):
            if r.url.path.endswith(resource + ".json"):
                data = json.loads(api.response(r).content)
                data[resource + "s"][0][resource + "_chembl_id"] = "CHEMBL999"
                return httpx.Response(200, json=data)

        api.hook = hook

    outcome, _, _ = chembl_case(configure)
    assert outcome.result.status == "failed"


def test_identical_duplicate_does_not_hide_missing_activity(chembl_case):
    def configure(api):
        def hook(r):
            if r.url.path.endswith("activity.json") and r.url.params["offset"] == "2":
                data = json.loads(api.response(r).content)
                data["activities"][0] = copy.deepcopy(api.data["activities"][0])
                return httpx.Response(200, json=data)

        api.hook = hook

    outcome, _, _ = chembl_case(configure)
    assert outcome.result.error.kind == "schema"
    assert "incomplete" in outcome.result.error.message
    assert outcome.result.counts["candidate_activities"] == 4


def test_service_failure_after_valid_page_retains_raw(chembl_case, tmp_path):
    def configure(api):
        api.hook = lambda r: (
            httpx.Response(503, content=b"down")
            if r.url.path.endswith("activity.json") and r.url.params["offset"] == "2"
            else None
        )

    outcome, _, _ = chembl_case(configure)
    assert outcome.result.error.kind == "http" and outcome.result.records == []
    assert outcome.result.counts["candidate_activities"] == 2
    folder = save_results(outcome, root=tmp_path)
    manifest = json.loads((folder / "sources/chembl-requests.json").read_text())
    assert len([e for e in manifest["exchanges"] if e["http_status"] == 503]) == 3
    assert any(
        (folder / "sources" / e["raw_path"]).read_bytes() == b"down"
        for e in manifest["exchanges"]
        if e["raw_path"]
    )


def test_last_transport_failure_does_not_reuse_previous_body(chembl_case):
    def configure(api):
        requests = []

        def hook(r):
            if r.url.path.endswith("target.json"):
                requests.append(1)
                if len(requests) == 1:
                    return httpx.Response(503, content=b"previous body")
                raise httpx.ReadTimeout("slow")

        api.hook = hook

    outcome, _, _ = chembl_case(configure)
    assert outcome.result.error.kind == "transport"
    assert outcome.result.error.http_status is None
    final = [e for e in outcome.exchanges if e.role == "targets"][-1]
    assert final.body is None and final.http_status is None


def test_protocol_error_does_not_retry(chembl_case):
    def configure(api):
        def hook(r):
            if r.url.path.endswith("target.json"):
                raise httpx.RemoteProtocolError("broken")

        api.hook = hook

    outcome, api, _ = chembl_case(configure)
    assert outcome.result.error.kind == "transport"
    assert len([r for r in api.requests if r.url.path.endswith("target.json")]) == 1


def test_non_up_status_and_unknown_version_warnings_are_ordered(chembl_case):
    def configure(api):
        api.hook = lambda r: (
            httpx.Response(200, json={"status": "DOWN"})
            if r.url.path.endswith("status.json")
            else None
        )

    outcome, _, _ = chembl_case(configure)
    assert outcome.result.status == "success"
    assert [(w["role"], w["code"]) for w in outcome.result.warnings] == [
        ("status_before", "release_metadata_unavailable"),
        ("status_before", "release_service_not_up"),
        ("status_after", "release_metadata_unavailable"),
        ("status_after", "release_service_not_up"),
    ]


def test_interrupt_propagates_without_final_probe(chembl_case):
    calls = []

    def configure(api):
        def hook(r):
            calls.append(str(r.url))
            if r.url.path.endswith("target.json"):
                raise KeyboardInterrupt

        api.hook = hook

    with pytest.raises(KeyboardInterrupt):
        chembl_case(configure)
    assert len(calls) == 2


@pytest.mark.parametrize(
    "bad", ["counts", "type", "values", "audit", "warning", "policy"]
)
def test_v2_reader_rejects_malformed_fields(chembl_case, bad):
    outcome, _, _ = chembl_case()
    data = asdict(outcome.result)
    if bad == "counts":
        data["counts"]["unique_ligands"] = 999
    elif bad == "type":
        data["records"][0]["target_confidence"] = True
    elif bad == "values":
        data["records"][0]["source_values"]["standard"]["value"] = 12
    elif bad == "audit":
        data["selection_audit"]["activities"][0]["included"] = "true"
    elif bad == "warning":
        data["warnings"] = [{"code": 1, "message": "bad", "role": None}]
    else:
        data["policy"] = "unknown"
    with pytest.raises(ValueError):
        source_result_from_dict(data)


@pytest.mark.parametrize("failure_mode", ["none", "bindingdb", "chembl", "both"])
def test_multi_source_cli_success_and_independent_failures(
    chembl_data, now, payload, tmp_path, monkeypatch, failure_mode
):
    from ligandscout import cli

    clock = Clock()
    api = API(copy.deepcopy(chembl_data), clock)
    if failure_mode in ("chembl", "both"):
        api.hook = lambda r: (
            httpx.Response(404) if r.url.path.endswith("target.json") else None
        )
    real_client = httpx.Client
    requests = []

    def handler(r):
        requests.append(r)
        if r.url.host == "www.ebi.ac.uk":
            return api(r)
        if failure_mode in ("bindingdb", "both"):
            return httpx.Response(404)
        return httpx.Response(200, json=payload)

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        cli.httpx,
        "Client",
        lambda **kwargs: real_client(transport=httpx.MockTransport(handler)),
    )
    monkeypatch.setattr(cli.time, "sleep", clock.sleep)
    monkeypatch.setattr(cli.time, "monotonic", clock.monotonic)
    response = CliRunner().invoke(cli.app, ["P00533"])
    assert response.exit_code == (0 if failure_mode == "none" else 1), response.output
    assert any(r.url.host == "www.ebi.ac.uk" for r in requests)
    report = (tmp_path / "results/P00533/report.md").read_text()
    assert ("| BindingDB | Unavailable |" in report) == (
        failure_mode in ("bindingdb", "both")
    )
    assert ("| ChEMBL | Unavailable |" in report) == (
        failure_mode in ("chembl", "both")
    )
    if failure_mode not in ("chembl", "both"):
        assert "| ChEMBL | 2 |" in report
    assert "chembl-ligands.html" in response.stdout


@pytest.mark.parametrize(
    "fail_path",
    ["chembl-requests.json", "chembl-ligands.html", "chembl.json", "report.md"],
)
def test_multi_source_publish_failure_cleans_temps(
    chembl_case, tmp_path, monkeypatch, fail_path
):
    from ligandscout import storage

    outcome, _, _ = chembl_case()
    real_replace = storage.os.replace

    def fail(source, destination):
        if destination.name == fail_path:
            raise OSError("disk full")
        return real_replace(source, destination)

    monkeypatch.setattr(storage.os, "replace", fail)
    with pytest.raises(OSError):
        save_results(outcome, root=tmp_path)
    assert not list(tmp_path.rglob(".ligandscout-*"))


def test_multi_source_staging_failure_and_interrupt_preserve_old_files(
    chembl_case, tmp_path, monkeypatch
):
    from ligandscout import storage

    outcome, _, _ = chembl_case()
    folder = save_results(outcome, root=tmp_path)
    before = {p: p.read_bytes() for p in folder.rglob("*") if p.is_file()}
    real_temp = storage.tempfile.NamedTemporaryFile
    calls = []

    def fail(**kwargs):
        calls.append(1)
        if len(calls) == 3:
            raise KeyboardInterrupt
        return real_temp(**kwargs)

    monkeypatch.setattr(storage.tempfile, "NamedTemporaryFile", fail)
    with pytest.raises(KeyboardInterrupt):
        save_results(outcome, root=tmp_path)
    assert all(p.read_bytes() == body for p, body in before.items())
    assert not list(tmp_path.rglob(".ligandscout-*"))


@pytest.mark.parametrize("activity_document", ["CHEMBL600", None])
def test_activity_document_precedence_and_assay_fallback(
    chembl_case, activity_document
):
    def configure(api):
        api.data["activities"][0]["document_chembl_id"] = activity_document
        api.data["assays"][0]["document_chembl_id"] = "CHEMBL999"
        api.data["documents"].append(
            {
                "document_chembl_id": "CHEMBL999",
                "pubmed_id": 999,
                "doi": "10.1000/assay-protocol",
            }
        )

    outcome, api, _ = chembl_case(configure)
    assert outcome.result.status == "success", outcome.result.error
    row = next(r for r in outcome.result.records if r.source_record_id == "1")
    expected = activity_document or "CHEMBL999"
    assert row.document_id == expected
    assert row.pmid == ("12345" if activity_document else "999")
    assert row.doi == (
        "10.1000/example" if activity_document else "10.1000/assay-protocol"
    )
    requested = {
        identifier
        for r in api.requests
        if r.url.path.endswith("document.json")
        for identifier in r.url.params["document_chembl_id__in"].split(",")
    }
    assert ("CHEMBL999" in requested) == (activity_document is None)


def test_real_abl1_activity_and_assay_can_have_distinct_documents(
    chembl_case, chembl_data, tmp_path
):
    pair = json.loads(
        (Path(__file__).parent / "fixtures/chembl-distinct-documents.json").read_text()
    )
    activity, assay = pair["activity"], pair["assay"]
    data = copy.deepcopy(chembl_data)
    data["targets"][0].update(
        target_chembl_id=activity["target_chembl_id"], pref_name="ABL1"
    )
    data["targets"][0]["target_components"][0]["accession"] = "P00519"
    data["targets"] = data["targets"][:1]
    data["activities"], data["assays"] = [activity], [assay]
    data["molecules"] = [{"molecule_chembl_id": activity["molecule_chembl_id"]}]
    data["documents"] = [
        {
            "document_chembl_id": activity["document_chembl_id"],
            "pubmed_id": 12345,
            "doi": "10.1000/activity-reference",
        }
    ]
    outcome, _, _ = chembl_case(data=data, accession="P00519")
    assert outcome.result.status == "success", outcome.result.error
    row = outcome.result.records[0]
    assert row.source_record_id == "23317250"
    assert row.target_accessions == ["P00519"]
    assert row.assay_id == "CHEMBL4879556"
    assert row.document_id == "CHEMBL4507256"
    assert row.variant_mutation == "E255K"
    assert row.pmid == "12345"
    assert outcome.result.counts["unique_ligands"] == 1
    page = render_html(outcome.result)
    assert "/document/CHEMBL4507256" in page
    assert "/document/CHEMBL6193874" not in page
    folder = save_results(outcome, root=tmp_path)
    manifest = json.loads((folder / "sources/chembl-requests.json").read_text())
    raw_assay = next(e for e in manifest["exchanges"] if e["role"] == "assays")
    raw = json.loads((folder / "sources" / raw_assay["raw_path"]).read_bytes())
    assert raw["assays"][0]["document_chembl_id"] == "CHEMBL6193874"


def test_different_but_invalid_assay_document_still_fails(chembl_case):
    def configure(api):
        api.data["assays"][0]["document_chembl_id"] = "invalid"

    outcome, _, _ = chembl_case(configure)
    assert outcome.result.error.kind == "schema"
    assert outcome.result.records == []
