"""UniProt annotation semantics, independent sources, and saved v4 results."""

import copy
import csv
import hashlib
import io
import json
from dataclasses import asdict

import httpx
import pytest

from ligandscout.api_http import Audit, request_json
from ligandscout.models import source_result_from_dict
from ligandscout.report import (
    render_html,
    render_index,
    render_report,
    render_summary,
)
from ligandscout.search import run_simple
from ligandscout.storage import save_results
from ligandscout.uniprot import fetch_uniprot


@pytest.fixture
def fetch_uni(now, uniprot_entry):
    def run(entry=None, handler=None, query="P00519"):
        sleeps = []

        def default(request):
            if request.url.host == "rest.uniprot.org":
                return httpx.Response(
                    200,
                    json=entry if entry is not None else uniprot_entry,
                    headers={"X-UniProt-Release": "2026_04", "Set-Cookie": "private"},
                )
            return httpx.Response(204)

        with httpx.Client(transport=httpx.MockTransport(handler or default)) as client:
            result = fetch_uniprot(
                query, client=client, sleep=sleeps.append, now=now, monotonic=lambda: 0
            )
            assert not client.is_closed
        return result, sleeps

    return run


def test_atp_cofactor_sites_and_reaction_not_ligands(fetch_uni):
    outcome, _ = fetch_uni()
    result = outcome.result
    assert result.status == "success" and result.source_version == "2026_04"
    assert [r.compound_id for r in result.records] == ["CHEBI:18420", "CHEBI:30616"]
    mg, atp = result.records
    assert len(atp.binding_sites) == 3
    assert atp.binding_sites[1]["residue_sequence"] == "K"
    assert atp.binding_sites[0]["location"]["start"]["value"] == 248
    assert atp.annotation_evidence[0]["category"] == "experimental"
    assert mg.ligand_roles == ["cofactor"] and mg.binding_sites == []
    assert mg.annotation_evidence[0]["category"] == "similarity"
    assert atp.measurement_value is None and atp.pmid is None
    assert atp.structure_provenance["status"] == "no_exact_mapping"
    assert atp.smiles is None and "CHEBI:999" not in [
        r.compound_id for r in result.records
    ]
    assert result.counts["unique_ligands"] == 2
    assert result.counts["binding_site_annotations"] == 3
    assert source_result_from_dict(asdict(result)) == result


def test_group_same_id_across_cofactor_and_site(fetch_uni, uniprot_entry):
    uniprot_entry["comments"][0]["cofactors"].append(
        {
            "name": "ATP renamed",
            "cofactorCrossReference": {"database": "ChEBI", "id": "CHEBI:30616"},
        }
    )
    outcome, _ = fetch_uni(uniprot_entry)
    atp = outcome.result.records[1]
    assert atp.ligand_roles == ["binding_site", "cofactor"]
    assert atp.source_values["aliases"] == ["ATP", "ATP renamed"]
    assert outcome.result.counts["cofactor_annotations"] == 2


def test_generic_context_is_not_a_ligand(fetch_uni, uniprot_entry):
    uniprot_entry["features"] = [uniprot_entry["features"][0]]
    uniprot_entry["features"][0]["ligand"] = {"name": "substrate"}
    uniprot_entry["comments"] = [uniprot_entry["comments"][1]]
    outcome, _ = fetch_uni(uniprot_entry)
    assert outcome.result.status == "no_matches" and outcome.result.records == []
    assert len(outcome.result.source_context["generic_sites"]) == 1
    page = render_html(outcome.result)
    assert "Generic binding sites" in page and "Catalytic Activity" in page
    assert "Reaction participants are context" in page
    assert "No specifically identified ligands" in page
    assert source_result_from_dict(asdict(outcome.result)) == outcome.result


def test_name_only_hash_and_id_differences(fetch_uni, uniprot_entry):
    first = uniprot_entry["features"][0]
    first["ligand"] = {"name": "  Test  molecule "}
    second = copy.deepcopy(first)
    second["ligand"] = {"name": "test molecule"}
    third = copy.deepcopy(first)
    third["ligand"] = {"name": "test molecule", "id": "CHEBI:123"}
    uniprot_entry.update(features=[first, second, third], comments=[])
    outcome, _ = fetch_uni(uniprot_entry)
    assert len(outcome.result.records) == 2
    record = outcome.result.records[1]
    assert (
        record.compound_id
        == "UNIPROT:P00519:NAME:" + hashlib.sha256(b"test molecule").hexdigest()
    )
    assert len(record.binding_sites) == 2
    assert record.structure_provenance["status"] == "not_applicable"
    assert len(outcome.exchanges) == 2  # core and one chemical search


def test_heme_ligand_part_is_not_a_third_ligand(fetch_uni, uniprot_entry):
    feature = uniprot_entry["features"][0]
    feature["ligand"] = {"name": "heme b", "id": "ChEBI:CHEBI:60344"}
    feature["ligandPart"] = {"name": "Fe", "id": "ChEBI:CHEBI:18248"}
    feature["featureCrossReferences"] = [
        {"database": "ChEBI", "id": "CHEBI:60344"},
        {"database": "ChEBI", "id": "CHEBI:18248"},
    ]
    oxygen = copy.deepcopy(feature)
    oxygen.pop("ligandPart")
    oxygen["ligand"] = {"name": "O2", "id": "CHEBI:15379"}
    oxygen.pop("featureCrossReferences")
    uniprot_entry.update(features=[feature, oxygen], comments=[])
    outcome, _ = fetch_uni(uniprot_entry)
    assert {r.compound_id for r in outcome.result.records} == {
        "CHEBI:60344",
        "CHEBI:15379",
    }
    heme = outcome.result.records[1]
    assert heme.binding_sites[0]["ligand_part"]["id"] == "CHEBI:18248"


@pytest.mark.parametrize(
    "sequence,modifier,expected",
    [(None, "UNSURE", None), ("P00519-2", "EXACT", None), (None, "EXACT", "K")],
)
def test_uncertain_and_alternate_sites(
    fetch_uni, uniprot_entry, sequence, modifier, expected
):
    f = uniprot_entry["features"][1]
    f["location"]["start"]["modifier"] = modifier
    f["location"]["end"]["modifier"] = modifier
    if sequence:
        f["location"]["sequence"] = sequence
    uniprot_entry.update(features=[f], comments=[])
    outcome, _ = fetch_uni(uniprot_entry)
    assert outcome.result.records[0].binding_sites[0]["residue_sequence"] == expected
    assert source_result_from_dict(asdict(outcome.result)) == outcome.result
    if sequence:
        assert "alternate sequence; residue letters unavailable" in render_html(
            outcome.result
        )


def test_polyprotein_containment(fetch_uni, uniprot_entry):
    def product(start, end, name):
        return {
            "type": "Chain",
            "description": name,
            "location": {"start": {"value": start}, "end": {"value": end}},
        }

    uniprot_entry["features"] += [
        product(200, 300, "Product A"),
        product(250, 350, "Product B"),
    ]
    outcome, _ = fetch_uni(uniprot_entry)
    sites = outcome.result.records[1].binding_sites
    assert [p["name"] for p in sites[0]["processed_products"]] == ["Product A"]
    assert [p["name"] for p in sites[1]["processed_products"]] == [
        "Product A",
        "Product B",
    ]


@pytest.mark.parametrize(
    "mutation",
    [
        "identity",
        "unreviewed",
        "sequence",
        "organism",
        "coordinate",
        "evidence",
        "type",
        "cofactor",
    ],
)
def test_core_validation_failures(fetch_uni, uniprot_entry, mutation):
    if mutation == "identity":
        uniprot_entry["primaryAccession"] = "P00699"
    if mutation == "unreviewed":
        uniprot_entry["entryType"] = "UniProtKB unreviewed (TrEMBL)"
    if mutation == "sequence":
        uniprot_entry["sequence"]["length"] = True
    if mutation == "organism":
        uniprot_entry["organism"]["taxonId"] = False
    if mutation == "coordinate":
        uniprot_entry["features"][0]["location"]["end"]["value"] = 999
    if mutation == "evidence":
        uniprot_entry["features"][0]["evidences"][0]["evidenceCode"] = "not evidence"
    if mutation == "type":
        uniprot_entry["features"][0]["type"] = None
    if mutation == "cofactor":
        uniprot_entry["comments"][0]["cofactors"] = None
    outcome, _ = fetch_uni(uniprot_entry)
    assert outcome.result.status == "failed" and outcome.result.records == []
    assert outcome.result.error.kind == "schema" and outcome.result.source_context == {}
    assert source_result_from_dict(asdict(outcome.result)) == outcome.result


def test_secondary_accession_resolution(fetch_uni, uniprot_entry):
    uniprot_entry["secondaryAccessions"] = ["P00699"]
    outcome, _ = fetch_uni(uniprot_entry, query="P00699")
    assert outcome.result.accession == "P00699"
    assert outcome.result.records[0].target_id == "P00519"
    assert "Resolved primary accession: P00519" in render_html(outcome.result)


def test_structural_failure_does_not_fail_core(fetch_uni, uniprot_entry):
    def handler(request):
        if request.url.host == "rest.uniprot.org":
            return httpx.Response(200, json=uniprot_entry)
        return httpx.Response(503)

    outcome, _ = fetch_uni(handler=handler)
    assert outcome.result.status == "success"
    assert all(
        r.structure_provenance["status"] == "unavailable"
        for r in outcome.result.records
    )
    assert outcome.result.warnings


def test_audited_storage_and_report_links(fetch_uni, tmp_path):
    outcome, _ = fetch_uni()
    folder = save_results(outcome, root=tmp_path)
    manifest = json.loads((folder / "sources/uniprot-requests.json").read_text())
    assert manifest["schema_version"] == 1 and len(manifest["exchanges"]) == 3
    for entry in manifest["exchanges"]:
        body = (folder / "sources" / entry["raw_path"]).read_bytes()
        assert hashlib.sha256(body).hexdigest() == entry["sha256"]
        assert "set-cookie" not in entry["response_headers"]
    result = source_result_from_dict(
        json.loads((folder / "sources/uniprot.json").read_text())
    )
    assert result == outcome.result
    index = (folder / "index.html").read_text()
    assert 'href="uniprot-ligands.html"' in index and "Records</th>" in index
    assert "Measurements:" not in render_summary(result, folder=folder)
    assert "Ligand annotation records: 2" in render_summary(result, folder=folder)
    page = (folder / "uniprot-ligands.html").read_text()
    assert "271 (K)" in page and "Binding site not annotated" in page
    assert "Cofactor context" in page and "Required for activity." in page
    assert "ligands.tsv" in page and 'href="index.html"' in page
    rows = list(
        csv.DictReader(
            io.StringIO((folder / "ligands.tsv").read_text()), delimiter="\t"
        )
    )
    assert all(
        r["measurement_value"] == "" and r["measurement_type"] == "" for r in rows
    )
    assert json.loads(rows[1]["binding_sites"])[1]["residue_sequence"] == "K"
    assert "Annotation context" in (folder / "report.md").read_text()


def test_stale_audit_cleanup_is_scoped(fetch_uni, tmp_path):
    outcome, _ = fetch_uni()
    folder = save_results(outcome, root=tmp_path)
    raw = folder / "sources/uniprot-raw"
    (raw / "999999.json").write_text("stale")
    (raw / "notes.txt").write_text("user")
    outcome.exchanges = outcome.exchanges[:1]
    save_results(outcome, root=tmp_path)
    assert (raw / "notes.txt").read_text() == "user"
    assert not (raw / "999999.json").exists() and not (raw / "000002.txt").exists()


@pytest.mark.parametrize("bad", ["count", "site", "evidence", "target", "assay"])
def test_saved_annotations_reject_corruption(fetch_uni, bad):
    outcome, _ = fetch_uni()
    data = asdict(outcome.result)
    if bad == "count":
        data["counts"]["binding_site_annotations"] = False
    if bad == "site":
        data["records"][1]["binding_sites"][1]["residue_sequence"] = "A"
    if bad == "evidence":
        data["records"][1]["annotation_evidence"][0]["category"] = "similarity"
    if bad == "target":
        data["records"][1]["target_id"] = "P00699"
    if bad == "assay":
        data["records"][1]["measurement_value"] = "1"
    with pytest.raises(ValueError):
        source_result_from_dict(data)


def test_annotation_html_escapes_content(fetch_uni, uniprot_entry):
    uniprot_entry["features"][0]["ligand"]["name"] = "<script>bad</script>"
    outcome, _ = fetch_uni(uniprot_entry)
    page = render_html(outcome.result)
    assert "<script>" not in page and "&lt;script&gt;" in page


@pytest.mark.parametrize(
    "binding_fails,uniprot_fails",
    [(False, False), (True, False), (False, True), (True, True)],
)
def test_sources_are_independent(
    payload, empty_uniprot, now, binding_fails, uniprot_fails
):
    def handler(request):
        if request.url.host == "www.bindingdb.org":
            return (
                httpx.Response(404)
                if binding_fails
                else httpx.Response(200, json=payload)
            )
        assert request.url.host == "rest.uniprot.org"
        return (
            httpx.Response(404)
            if uniprot_fails
            else httpx.Response(200, json=empty_uniprot)
        )

    with httpx.Client(transport=httpx.MockTransport(handler)) as client:
        run = run_simple(
            "P00699", client=client, sleep=lambda _: None, now=now, monotonic=lambda: 0
        )
    assert [s.result.status for s in run.sources] == [
        "failed" if binding_fails else "success",
        "failed" if uniprot_fails else "no_matches",
    ]
    assert "BindingDB" in render_index(run) and "UniProt" in render_index(run)
    assert "## UniProt" in render_report(run)


def component(
    identifier="ATP",
    chemical="CHEBI:30616",
    smiles="C[C@H](N)O",
    field="SMILES_stereo",
    key="KEY",
):
    return {
        "rcsb_id": identifier,
        "chem_comp": {"id": identifier},
        "rcsb_chem_comp_related": [
            {"resource_name": "ChEBI", "resource_accession_code": chemical}
        ],
        "rcsb_chem_comp_descriptor": {field: smiles, "InChIKey": key},
    }


@pytest.mark.parametrize(
    "mode,expected",
    [
        ("stereo", "available"),
        ("plain", "available"),
        ("mismatch", "no_exact_mapping"),
        ("wrong_resource", "no_exact_mapping"),
        ("no_descriptor", "no_descriptor"),
        ("conflict", "ambiguous"),
        ("identical", "available"),
        ("broken", "unavailable"),
        ("wrong_identity", "unavailable"),
    ],
)
def test_ccd_identity_and_descriptor_selection(
    fetch_uni, uniprot_entry, mode, expected
):
    uniprot_entry.update(features=[uniprot_entry["features"][0]], comments=[])
    candidates = ["ATP", "OTHER"] if mode in ("conflict", "identical") else ["ATP"]

    def handler(request):
        if request.url.host == "rest.uniprot.org":
            return httpx.Response(200, json=uniprot_entry)
        if request.url.host == "search.rcsb.org":
            query = json.loads(request.url.params["json"])
            assert query["query"]["parameters"]["value"] == "CHEBI:30616"
            return httpx.Response(
                200,
                json={
                    "result_type": "mol_definition",
                    "total_count": len(candidates),
                    "result_set": [{"identifier": c} for c in candidates],
                },
            )
        identifier = request.url.path.rsplit("/", 1)[1]
        data = component(
            identifier=identifier,
            chemical="CHEBI:15422" if mode == "mismatch" else "CHEBI:30616",
            field="SMILES" if mode == "plain" else "SMILES_stereo",
        )
        if mode == "wrong_resource":
            data["rcsb_chem_comp_related"][0]["resource_name"] = "Other"
        if mode == "no_descriptor":
            data.pop("rcsb_chem_comp_descriptor")
        if mode == "conflict" and identifier == "OTHER":
            data["rcsb_chem_comp_descriptor"]["SMILES_stereo"] = "CC"
        if mode == "broken":
            data["rcsb_chem_comp_descriptor"] = []
        if mode == "wrong_identity":
            data["rcsb_id"] = "OTHER"
        return httpx.Response(200, json=data)

    outcome, _ = fetch_uni(handler=handler)
    record = outcome.result.records[0]
    assert (
        outcome.result.status == "success"
        and record.structure_provenance["status"] == expected
    )
    assert bool(record.smiles) == (expected == "available")
    assert source_result_from_dict(asdict(outcome.result)) == outcome.result
    page = render_html(outcome.result)
    if expected == "available":
        assert '<code class="smiles">' in page
        if mode == "plain":
            assert "non-stereospecific fallback" in page
    if mode == "mismatch":
        assert record.smiles is None


@pytest.mark.parametrize(
    "issue",
    ["duplicate", "changed_total", "short_page", "late_204", "late_failure", "valid"],
)
def test_ccd_pagination_and_partial_failures(fetch_uni, uniprot_entry, issue):
    uniprot_entry.update(features=[uniprot_entry["features"][0]], comments=[])

    def handler(request):
        if request.url.host == "rest.uniprot.org":
            return httpx.Response(200, json=uniprot_entry)
        if request.url.host == "search.rcsb.org":
            start = json.loads(request.url.params["json"])["request_options"][
                "paginate"
            ]["start"]
            if start and issue == "late_204":
                return httpx.Response(204)
            if start and issue == "late_failure":
                return httpx.Response(404)
            ids = [f"C{i}" for i in range(100)] if not start else ["C100"]
            if issue == "short_page" and not start:
                ids = ids[:5]
            if issue == "duplicate" and start:
                ids = ["C0"]
            total = 102 if issue == "changed_total" and start else 101
            return httpx.Response(
                200,
                json={
                    "result_type": "mol_definition",
                    "total_count": total,
                    "result_set": [{"identifier": i} for i in ids],
                },
            )
        return httpx.Response(200, json=component(request.url.path.rsplit("/", 1)[1]))

    outcome, _ = fetch_uni(handler=handler)
    record = outcome.result.records[0]
    assert record.structure_provenance["status"] == (
        "available" if issue == "valid" else "unavailable"
    )
    assert outcome.result.status == "success"
    if issue == "valid":
        assert len(record.structure_provenance["candidate_components"]) == 101


@pytest.mark.parametrize("body", [b"", b'""', b" \n"])
def test_empty_entry_retries_are_indeterminate(fetch_uni, body):
    calls = []

    def handler(request):
        calls.append(request)
        return httpx.Response(200, content=body)

    outcome, _ = fetch_uni(handler=handler)
    assert len(calls) == 3 and outcome.result.status == "failed"
    assert "indeterminate" in outcome.result.error.message
    assert [e.attempt for e in outcome.exchanges] == [1, 2, 3]


@pytest.mark.parametrize("status", [404, 204])
def test_core_http_failures_are_not_empty_annotations(fetch_uni, status):
    outcome, _ = fetch_uni(handler=lambda _: httpx.Response(status))
    assert len(outcome.exchanges) == 1 and outcome.result.status == "failed"
    assert outcome.result.error.http_status == status


def test_audit_transport_failures_and_retry_after(now):
    ticks = [0.0]

    def sleep(delay):
        ticks[0] += delay

    starts = []

    def handler(request):
        starts.append(ticks[0])
        if len(starts) == 1:
            raise httpx.ConnectError("no connection")
        if len(starts) == 2:
            return httpx.Response(429, headers={"Retry-After": "3"})
        return httpx.Response(200, json={"ok": True})

    audit = Audit()
    with httpx.Client(transport=httpx.MockTransport(handler)) as client:
        result = request_json(
            "https://rest.uniprot.org/test",
            params={},
            role="entry",
            chemical_id=None,
            audit=audit,
            client=client,
            sleep=sleep,
            now=now,
            monotonic=lambda: ticks[0],
        )
    assert result.error is None and starts == [0, 1, 4]
    assert [e.request_id for e in audit.exchanges] == [1, 1, 1]
    assert (
        audit.exchanges[0].body is None
        and audit.exchanges[0].transport_error == "ConnectError"
    )


def test_audit_redirect_rate_limit_and_policy(now):
    ticks = [0.0]
    starts = []

    def sleep(delay):
        ticks[0] += delay

    def handler(request):
        starts.append(ticks[0])
        return (
            httpx.Response(302, headers={"Location": "/new"})
            if len(starts) == 1
            else httpx.Response(302, headers={"Location": "http://unsafe.example/"})
        )

    audit = Audit()
    with httpx.Client(transport=httpx.MockTransport(handler)) as client:
        result = request_json(
            "https://rest.uniprot.org/test",
            params={},
            role="entry",
            chemical_id=None,
            audit=audit,
            client=client,
            sleep=sleep,
            now=now,
            monotonic=lambda: ticks[0],
        )
    assert result.error.kind == "http" and starts == [0, 0.5]
    assert [e.redirect_hop for e in audit.exchanges] == [0, 1]


def test_long_retry_after_stops(now):
    audit = Audit()
    sleeps = []
    with httpx.Client(
        transport=httpx.MockTransport(
            lambda _: httpx.Response(503, headers={"Retry-After": "61"})
        )
    ) as client:
        result = request_json(
            "https://rest.uniprot.org/test",
            params={},
            role="entry",
            chemical_id=None,
            audit=audit,
            client=client,
            sleep=sleeps.append,
            now=now,
            monotonic=lambda: 0,
        )
    assert not sleeps and len(audit.exchanges) == 1 and result.error


@pytest.mark.parametrize("bad", [None, 3, [], True])
def test_malformed_primary_is_a_source_failure(fetch_uni, uniprot_entry, bad):
    uniprot_entry["primaryAccession"] = bad
    outcome, _ = fetch_uni(uniprot_entry)
    assert outcome.result.status == "failed" and outcome.result.error.kind == "schema"


def test_unknown_positions_missing_evidence_and_unknown_features(
    fetch_uni, uniprot_entry
):
    feature = uniprot_entry["features"][0]
    feature["location"] = {}
    feature.pop("evidences")
    uniprot_entry.update(
        features=[feature, {"type": "Future annotation", "anything": 42}], comments=[]
    )
    outcome, _ = fetch_uni(uniprot_entry)
    binding = outcome.result.records[0].binding_sites[0]
    assert binding["location"]["start"] == {"value": None, "modifier": "UNKNOWN"}
    assert binding["residue_sequence"] is None and not binding["evidence"]
    assert 'href="sources/uniprot.json"' in render_html(outcome.result)


def test_optional_context_failure_and_malformed_cofactor_note(fetch_uni, uniprot_entry):
    uniprot_entry["comments"].append({"commentType": "FUNCTION", "texts": "bad"})
    outcome, _ = fetch_uni(uniprot_entry)
    assert outcome.result.status == "success" and any(
        w["code"] == "invalid_optional_context" for w in outcome.result.warnings
    )
    uniprot_entry["comments"][0]["note"] = {"texts": "bad"}
    outcome, _ = fetch_uni(uniprot_entry)
    assert outcome.result.status == "failed"


def test_request_redirect_keeps_location_query(now):
    seen = []

    def handler(request):
        seen.append(str(request.url))
        if len(seen) == 1:
            return httpx.Response(302, headers={"Location": "/new?token=abc"})
        assert request.url.params["token"] == "abc"
        return httpx.Response(200, json={"ok": True})

    audit = Audit()
    with httpx.Client(transport=httpx.MockTransport(handler)) as client:
        response = request_json(
            "https://rest.uniprot.org/test",
            params={},
            role="entry",
            chemical_id=None,
            audit=audit,
            client=client,
            sleep=lambda _: None,
            now=now,
            monotonic=lambda: 0,
        )
    assert response.error is None and audit.exchanges[-1].parameters == {"token": "abc"}


def test_partial_ccd_candidates_never_supply_selected_smiles(fetch_uni, uniprot_entry):
    uniprot_entry.update(features=[uniprot_entry["features"][0]], comments=[])

    def handler(request):
        if request.url.host == "rest.uniprot.org":
            return httpx.Response(200, json=uniprot_entry)
        if request.url.host == "search.rcsb.org":
            return httpx.Response(
                200,
                json={
                    "result_type": "mol_definition",
                    "total_count": 2,
                    "result_set": [{"identifier": "ATP"}, {"identifier": "OTHER"}],
                },
            )
        if request.url.path.endswith("ATP"):
            return httpx.Response(200, json=component())
        return httpx.Response(404)

    outcome, _ = fetch_uni(handler=handler)
    record = outcome.result.records[0]
    assert (
        record.smiles is None and record.structure_provenance["status"] == "unavailable"
    )
    assert len(record.structure_provenance["candidate_components"]) == 1


def test_empty_lookup_has_a_raw_audit_file(fetch_uni, tmp_path):
    outcome, _ = fetch_uni(handler=lambda _: httpx.Response(200, content=b""))
    folder = save_results(outcome, root=tmp_path)
    manifest = json.loads((folder / "sources/uniprot-requests.json").read_text())
    assert len(manifest["exchanges"]) == 3
    for entry in manifest["exchanges"]:
        assert (folder / "sources" / entry["raw_path"]).read_bytes() == b""
        assert entry["sha256"] == hashlib.sha256(b"").hexdigest()


def test_two_source_index_files_and_cli_failure_isolation(
    monkeypatch, tmp_path, payload, empty_uniprot
):
    from typer.testing import CliRunner

    from ligandscout import cli

    monkeypatch.chdir(tmp_path)
    original = httpx.Client

    def handler(request):
        if request.url.host == "www.bindingdb.org":
            return httpx.Response(404)
        assert request.url.host == "rest.uniprot.org"
        return httpx.Response(200, json=empty_uniprot)

    monkeypatch.setattr(
        cli.httpx,
        "Client",
        lambda **kwargs: original(transport=httpx.MockTransport(handler)),
    )
    result = CliRunner().invoke(cli.app, ["P00699"])
    assert result.exit_code == 1
    folder = tmp_path / "results/P00699"
    page = (folder / "index.html").read_text()
    for name in ("bindingdb", "uniprot"):
        assert (folder / f"{name}-ligands.html").is_file()
        assert (folder / f"sources/{name}.json").is_file()
        assert f'href="{name}-ligands.html"' in page
    assert "Source: BindingDB | Status: failed" in result.stdout
    assert "Source: UniProt | Status: no_matches" in result.stdout
    assert "Details: index.html" in result.stdout and "Ligands:" not in result.stdout


def test_optional_metadata_is_unknown_not_invented(fetch_uni, uniprot_entry):
    uniprot_entry["entryAudit"] = {
        "entryVersion": True,
        "sequenceVersion": -1,
        "lastAnnotationUpdateDate": "not a date",
    }
    uniprot_entry["references"] = {}
    outcome, _ = fetch_uni(uniprot_entry)
    entry = outcome.result.source_context["entry"]
    assert outcome.result.status == "success"
    assert entry["entry_audit"]["entryVersion"] is None
    assert entry["entry_audit"]["sequenceVersion"] is None
    assert entry["entry_audit"]["lastAnnotationUpdateDate"] is None
    assert entry["references"] == []
    assert (
        len(
            [
                w
                for w in outcome.result.warnings
                if w["code"] == "invalid_optional_metadata"
            ]
        )
        == 4
    )


def test_invalid_optional_reaction_is_not_mined_for_ligands(fetch_uni, uniprot_entry):
    uniprot_entry["comments"][1]["reaction"]["name"] = False
    outcome, _ = fetch_uni(uniprot_entry)
    assert len(outcome.result.records) == 2 and outcome.result.status == "success"
    assert all(
        c["comment_type"] != "CATALYTIC ACTIVITY"
        for c in outcome.result.source_context["context_comments"]
    )
    assert any(w["code"] == "invalid_optional_context" for w in outcome.result.warnings)


def test_context_summary_is_an_extract_preserving_qualifiers():
    from ligandscout.uniprot_report import summarize_text

    source = "Does not bind DNA (By similarity). May regulate growth. Interacts with ABL1. This fourth sentence must be omitted."
    summary, omitted = summarize_text([source])
    assert (
        summary
        == "Does not bind DNA (By similarity). May regulate growth. Interacts with ABL1."
    )
    assert omitted
    assert summarize_text(
        [
            "Kinase activity was measured e.g. By a fluorescence assay. A second sentence."
        ]
    )[0].startswith("Kinase activity was measured e.g. By a fluorescence assay.")


def test_html_binding_site_bullets_and_source_prose(fetch_uni, uniprot_entry):
    from html.parser import HTMLParser

    uniprot_entry["comments"].append(
        {
            "commentType": "FUNCTION",
            "texts": [
                {
                    "value": "May regulate growth (By similarity). Does not bind DNA. Interacts with ABL1. The fourth sentence is omitted."
                }
            ],
        }
    )
    outcome, _ = fetch_uni(uniprot_entry)
    page = render_html(outcome.result)
    assert "<th>Annotation details</th>" not in page
    assert "<pre>" not in page and "<summary>Annotation details</summary>" not in page
    assert "annotation_path" not in page and "candidate_components" not in page
    assert "May regulate growth (By similarity)." in page
    assert "The fourth sentence is omitted" not in page
    assert "First three source sentences" in page
    assert "https://www.uniprot.org/uniprotkb/P00519/entry" in page
    assert "https://rest.uniprot.org/uniprotkb/P00519.json" in page
    assert "Required for activity." in page

    class Cells(HTMLParser):
        def __init__(self):
            super().__init__()
            self.cells = []
            self.current = None

        def handle_starttag(self, tag, attrs):
            if tag == "td":
                self.current = []

        def handle_data(self, data):
            if self.current is not None:
                self.current.append(data)

        def handle_endtag(self, tag):
            if tag == "td" and self.current is not None:
                self.cells.append(" ".join(self.current))
                self.current = None

    parser = Cells()
    parser.feed(page)
    assert len(parser.cells) == 8  # four columns per ligand row
    assert "Binding site not annotated." in parser.cells[3]
    assert "271 (K)" in parser.cells[7]
    assert '<ul class="binding-sites"><li>248–256' in page
    assert "</li><li>271 (K)</li><li>316–322" in page
    assert "Cofactor annotations and notes" not in page


def test_context_formats_structured_fields_without_dumping_json(
    fetch_uni, uniprot_entry
):
    uniprot_entry["comments"].append(
        {
            "commentType": "BIOPHYSICOCHEMICAL PROPERTIES",
            "phDependence": {"texts": [{"value": "Optimum pH is 7."}]},
            "kineticParameters": {
                "michaelisConstants": [
                    {"constant": 2, "unit": "mM", "substrate": "ATP"}
                ]
            },
        }
    )
    outcome, _ = fetch_uni(uniprot_entry)
    page = render_html(outcome.result)
    assert "Optimum pH is 7." in page and "Michaelis constant: 2 mM for ATP." in page
    assert "phDependence" not in page and "<pre>" not in page
