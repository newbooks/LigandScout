"""Structural selection, optional enrichment and saved evidence regressions."""

import copy
import json
from dataclasses import asdict

import httpx
import pytest

from ligandscout.models import source_result_from_dict, upgrade_result
from ligandscout.rcsb import SEARCH, fetch_rcsb, mapping
from ligandscout.report import render_summary
from ligandscout.storage import save_results


@pytest.fixture
def rcsb_data():
    def contact(chain="A", pos=1, entity="1"):
        return {
            "target_asym_id": chain,
            "target_entity_id": entity,
            "target_seq_id": pos,
            "target_comp_id": "ALA",
            "target_auth_seq_id": -2,
            "distance": 3.2,
            "target_model_id": 1,
            "comp_id": "STI",
        }

    poly = {
        "rcsb_id": "1ABC_1",
        "rcsb_polymer_entity_container_identifiers": {
            "entry_id": "1ABC",
            "entity_id": "1",
            "asym_ids": ["A", "B"],
            "reference_sequence_identifiers": [
                {"database_name": "UniProt", "database_accession": "P00519"}
            ],
        },
        "entity_poly": {"rcsb_sample_sequence_length": 10},
        "rcsb_polymer_entity_align": [
            {
                "reference_database_name": "UniProt",
                "reference_database_accession": "P00519",
                "aligned_regions": [
                    {"entity_beg_seq_id": 1, "ref_beg_seq_id": 5, "length": 2}
                ],
            }
        ],
    }
    data = {
        "polymer_entity/1ABC/1": poly,
        "entry/1ABC": {
            "rcsb_id": "1ABC",
            "rcsb_entry_container_identifiers": {
                "entry_id": "1ABC",
                "non_polymer_entity_ids": ["2", "3", "4", "5"],
            },
            "rcsb_entry_info": {
                "structure_determination_methodology": "experimental",
                "resolution_combined": [2.1],
            },
            "exptl": [{"method": "X-RAY DIFFRACTION"}],
        },
    }
    for entity, chain, comp, rows in [
        ("2", "C", "STI", [contact("A"), contact("B")]),
        ("3", "D", "CL", [contact(pos=3) | {"comp_id": "CL"}]),
        ("4", "E", "HOH", []),
        ("5", "F", "XXX", [contact("Z", entity="99") | {"comp_id": "XXX"}]),
    ]:
        data[f"nonpolymer_entity/1ABC/{entity}"] = {
            "rcsb_id": f"1ABC_{entity}",
            "rcsb_nonpolymer_entity_container_identifiers": {
                "entry_id": "1ABC",
                "entity_id": entity,
                "asym_ids": [chain],
                "nonpolymer_comp_id": comp,
            },
            "pdbx_entity_nonpoly": {"comp_id": comp, "name": comp},
        }
        data[f"nonpolymer_entity_instance/1ABC/{chain}"] = {
            "rcsb_id": f"1ABC.{chain}",
            "rcsb_nonpolymer_entity_instance_container_identifiers": {
                "entry_id": "1ABC",
                "entity_id": entity,
                "asym_id": chain,
                "comp_id": comp,
                "auth_seq_id": "4",
                "auth_asym_id": "A",
            },
            "rcsb_target_neighbors": rows,
        }
    for comp in ("STI", "CL"):
        data["chemcomp/" + comp] = {
            "rcsb_id": comp,
            "chem_comp": {"id": comp, "name": comp},
            "rcsb_chem_comp_descriptor": {"SMILES_stereo": "CN", "InChIKey": "key"}
            if comp == "STI"
            else {},
        }
    return data


@pytest.fixture
def fetch_rcsb_fixture(now, rcsb_data):
    def call(data=None, *, handler=None, context=None):
        requests, sleeps = [], []
        objects = rcsb_data if data is None else data

        def default(request):
            requests.append(request)
            if str(request.url).startswith(SEARCH):
                q = json.loads(request.url.params["json"])
                assert q["request_options"]["results_content_type"] == ["experimental"]
                return httpx.Response(
                    200,
                    json={
                        "result_type": "polymer_entity",
                        "total_count": 1,
                        "result_set": [{"identifier": "1ABC_1"}],
                    },
                )
            key = request.url.path.removeprefix("/rest/v1/core/")
            value = objects[key]
            return (
                value
                if isinstance(value, httpx.Response)
                else httpx.Response(200, json=value)
            )

        with httpx.Client(transport=httpx.MockTransport(handler or default)) as client:
            outcome = fetch_rcsb(
                "P00519",
                client=client,
                sleep=sleeps.append,
                now=now,
                monotonic=lambda: 0,
                uniprot_context=context,
            )
            assert not client.is_closed
        return outcome, requests, sleeps

    return call


def test_baseline_roundtrip_and_reports(fetch_rcsb_fixture, tmp_path):
    out, requests, _ = fetch_rcsb_fixture()
    r = out.result
    assert r.status == "success"
    assert {
        k: r.counts[k]
        for k in (
            "candidate_entries",
            "accepted_target_entities",
            "examined_ligand_instances",
            "included_ligand_instances",
            "excluded_water_instances",
            "instances_without_reported_contacts",
            "instances_without_target_contacts",
            "unique_ligands",
            "contact_sites",
            "contact_residues",
            "mapped_contact_residues",
            "smiles_available",
            "smiles_unavailable",
        )
    } == {
        "candidate_entries": 1,
        "accepted_target_entities": 1,
        "examined_ligand_instances": 4,
        "included_ligand_instances": 2,
        "excluded_water_instances": 1,
        "instances_without_reported_contacts": 0,
        "instances_without_target_contacts": 1,
        "unique_ligands": 2,
        "contact_sites": 3,
        "contact_residues": 3,
        "mapped_contact_residues": 2,
        "smiles_available": 1,
        "smiles_unavailable": 1,
    }
    assert source_result_from_dict(asdict(r)) == r
    assert upgrade_result(r) is r
    folder = save_results(out, root=tmp_path)
    page = (folder / "rcsb-ligands.html").read_text()
    assert "P00519:5" in page and "outside_alignment" in page
    assert '<ul class="binding-sites">' in page and "<pre>" not in page
    assert "rcsb-ligands.html" in (folder / "index.html").read_text()
    assert "structure_associated_ligand" in (folder / "ligands.tsv").read_text()
    assert "Contact sites: 3" in render_summary(out, folder=folder)
    assert "Measurements:" not in render_summary(out, folder=folder)
    manifest = json.loads((folder / "sources/rcsb-requests.json").read_text())
    assert len(manifest["exchanges"]) == len(requests)
    for exchange in manifest["exchanges"]:
        assert (folder / "sources" / exchange["raw_path"]).exists()


@pytest.mark.parametrize(
    "issue",
    [
        "other_accession",
        "absent_refs",
        "wrong_database",
        "different_method",
        "absent_contacts",
        "unrelated_contacts",
    ],
)
def test_exclusions(fetch_rcsb_fixture, rcsb_data, issue):
    refs = rcsb_data["polymer_entity/1ABC/1"][
        "rcsb_polymer_entity_container_identifiers"
    ]
    if issue == "other_accession":
        refs["reference_sequence_identifiers"][0]["database_accession"] = "P00520"
    if issue == "absent_refs":
        refs.pop("reference_sequence_identifiers")
    if issue == "wrong_database":
        refs["reference_sequence_identifiers"][0]["database_name"] = "Other"
    if issue == "different_method":
        rcsb_data["entry/1ABC"]["rcsb_entry_info"][
            "structure_determination_methodology"
        ] = "computational"
    if issue in ("absent_contacts", "unrelated_contacts"):
        for chain in ("C", "D"):
            obj = rcsb_data["nonpolymer_entity_instance/1ABC/" + chain]
            if issue == "absent_contacts":
                obj.pop("rcsb_target_neighbors")
            else:
                for row in obj["rcsb_target_neighbors"]:
                    row.update(target_entity_id="99", target_asym_id="Z")
    out, requests, _ = fetch_rcsb_fixture()
    assert out.result.status == "no_matches"
    assert out.result.counts["unique_ligands"] == 0
    if issue in ("other_accession", "absent_refs", "wrong_database"):
        assert len(requests) == 2


@pytest.mark.parametrize(
    "issue",
    [
        "wrong_entry",
        "wrong_ligand",
        "wrong_parent",
        "null_contacts",
        "bad_contact",
        "bad_distance",
        "boolean_position",
        "mismatched_chain",
        "bad_entity_length",
    ],
)
def test_required_failure_drops_partial_records(fetch_rcsb_fixture, rcsb_data, issue):
    obj = rcsb_data["nonpolymer_entity_instance/1ABC/D"]
    row = obj["rcsb_target_neighbors"][0]
    if issue == "wrong_entry":
        rcsb_data["entry/1ABC"]["rcsb_id"] = "XXXX"
    if issue == "wrong_ligand":
        obj["rcsb_nonpolymer_entity_instance_container_identifiers"]["comp_id"] = "STI"
    if issue == "wrong_parent":
        obj["rcsb_nonpolymer_entity_instance_container_identifiers"]["entity_id"] = "99"
    if issue == "null_contacts":
        obj["rcsb_target_neighbors"] = None
    if issue == "bad_contact":
        row.pop("target_comp_id")
    if issue == "bad_distance":
        row["distance"] = -1
    if issue == "boolean_position":
        row["target_seq_id"] = True
    if issue == "mismatched_chain":
        row["target_entity_id"] = "99"
    if issue == "bad_entity_length":
        rcsb_data["polymer_entity/1ABC/1"]["entity_poly"][
            "rcsb_sample_sequence_length"
        ] = False
    out, _, _ = fetch_rcsb_fixture()
    assert out.result.status == "failed" and out.result.records == []
    assert out.result.counts["unique_ligands"] is None
    assert source_result_from_dict(asdict(out.result)) == out.result


@pytest.mark.parametrize(
    "issue",
    ["outage", "identity", "malformed", "no_descriptor", "fallback", "key_only"],
)
def test_optional_chemical_enrichment(fetch_rcsb_fixture, rcsb_data, issue):
    obj = rcsb_data["chemcomp/STI"]
    if issue == "outage":
        rcsb_data["chemcomp/STI"] = httpx.Response(503)
    if issue == "identity":
        obj["chem_comp"]["id"] = "OTHER"
    if issue == "malformed":
        obj["rcsb_chem_comp_descriptor"] = []
    if issue == "no_descriptor":
        obj.pop("rcsb_chem_comp_descriptor")
    if issue == "fallback":
        obj["rcsb_chem_comp_descriptor"] = {"SMILES": "CN"}
    if issue == "key_only":
        obj["rcsb_chem_comp_descriptor"] = {"InChIKey": "key"}
    out, _, _ = fetch_rcsb_fixture()
    assert out.result.status == "success"
    record = next(r for r in out.result.records if r.compound_id == "CCD:STI")
    assert record.smiles == ("CN" if issue == "fallback" else None)
    if issue == "key_only":
        assert record.inchikey == "key"


@pytest.mark.parametrize(
    "issue,reason",
    [
        ("missing", "no_alignment"),
        ("gap", "outside_alignment"),
        ("conflict", "conflicting_alignment"),
        ("malformed", "invalid_alignment"),
        ("bounds", "uniprot_out_of_bounds"),
    ],
)
def test_mapping_edge_cases(rcsb_data, issue, reason):
    entity = rcsb_data["polymer_entity/1ABC/1"]
    align = entity["rcsb_polymer_entity_align"]
    if issue == "missing":
        entity.pop("rcsb_polymer_entity_align")
    if issue == "conflict":
        align[0]["aligned_regions"].append(
            {"entity_beg_seq_id": 1, "ref_beg_seq_id": 6, "length": 2}
        )
    if issue == "malformed":
        align[0]["aligned_regions"][0]["length"] = True
    context = (
        {"primary_accession": "P00519", "sequence": "AAA"}
        if issue == "bounds"
        else None
    )
    result = mapping(entity, 3 if issue == "gap" else 1, ["P00519"], "P00519", context)
    assert result["mapping_reason"] == reason and result["uniprot_position"] is None


def test_sequence_mismatch_does_not_drop_position(fetch_rcsb_fixture):
    out, _, _ = fetch_rcsb_fixture(
        context={
            "primary_accession": "P00519",
            "secondary_accessions": [],
            "sequence": "GGGGGGGGGG",
        }
    )
    residue = out.result.records[1].binding_sites[0]["residues"][0]
    assert residue["uniprot_position"] == 5 and residue["sequence_check"] == "mismatch"


@pytest.mark.parametrize(
    "response",
    [
        httpx.Response(204),
        httpx.Response(200, json={"result_type": "polymer_entity", "total_count": 0}),
    ],
)
def test_empty_search(fetch_rcsb_fixture, response):
    out, _, _ = fetch_rcsb_fixture(handler=lambda _: response)
    assert (
        out.result.status == "no_matches"
        and out.result.selection_audit["outcome_reason"] == "no_exact_structures"
    )


@pytest.mark.parametrize(
    "issue", ["empty", "malformed", "duplicate", "short", "wrong_type", "http"]
)
def test_search_failures(fetch_rcsb_fixture, issue):
    body = {
        "result_type": "polymer_entity",
        "total_count": 1,
        "result_set": [{"identifier": "1ABC_1"}],
    }
    if issue == "duplicate":
        body.update(total_count=2, result_set=[{"identifier": "1ABC_1"}] * 2)
    if issue == "short":
        body["total_count"] = 200
    if issue == "wrong_type":
        body["result_type"] = "entry"
    response = (
        httpx.Response(200, content=b"")
        if issue == "empty"
        else httpx.Response(200, content=b"{")
        if issue == "malformed"
        else httpx.Response(404)
        if issue == "http"
        else httpx.Response(200, json=body)
    )
    out, _, sleeps = fetch_rcsb_fixture(handler=lambda _: response)
    assert out.result.status == "failed"
    if issue == "empty":
        assert len(out.exchanges) == 3 and 1 in sleeps and 2 in sleeps


def test_overwrite_and_html_escaping(fetch_rcsb_fixture, rcsb_data, tmp_path):
    rcsb_data["chemcomp/STI"]["chem_comp"]["name"] = "<script>alert(1)</script>"
    out, _, _ = fetch_rcsb_fixture()
    folder = save_results(out, root=tmp_path)
    assert "<script>" not in (folder / "rcsb-ligands.html").read_text()
    assert "&lt;script&gt;" in (folder / "rcsb-ligands.html").read_text()
    unrelated = folder / "sources/unrelated.txt"
    unrelated.write_text("keep")
    empty, _, _ = fetch_rcsb_fixture(handler=lambda _: httpx.Response(204))
    save_results(empty, root=tmp_path)
    assert "CCD:STI" not in (folder / "ligands.tsv").read_text()
    assert len(list((folder / "sources/rcsb-raw").iterdir())) == 1
    assert unrelated.read_text() == "keep"


@pytest.mark.parametrize(
    "issue", [None, "changed_total", "duplicate", "late_empty", "no_progress"]
)
def test_pagination_is_complete_and_audited(now, issue):
    from ligandscout.rcsb import Adapter, RCSBError

    offsets = []

    def handler(request):
        start = json.loads(request.url.params["json"])["request_options"]["paginate"][
            "start"
        ]
        offsets.append(start)
        if start and issue == "late_empty":
            return httpx.Response(204)
        total = 102 if start and issue == "changed_total" else 101
        rows = (
            [{"identifier": f"1ABC_{i}"} for i in range(1, 101)]
            if start == 0
            else [{"identifier": "1ABC_1" if issue == "duplicate" else "1ABC_101"}]
        )
        if start and issue == "no_progress":
            rows = []
        return httpx.Response(
            200,
            json={
                "result_type": "polymer_entity",
                "total_count": total,
                "result_set": rows,
            },
        )

    with httpx.Client(transport=httpx.MockTransport(handler)) as client:
        adapter = Adapter("P00519", client, lambda _: None, now, lambda: 0, None)
        if issue:
            with pytest.raises(RCSBError):
                adapter.search()
        else:
            assert len(adapter.search()) == 101
        assert offsets == [0, 100] and len(adapter.audit.exchanges) == 2


def test_secondary_resolution_searches_both_and_deduplicates(fetch_rcsb_fixture):
    context = {
        "primary_accession": "P00520",
        "secondary_accessions": ["P00519"],
        "sequence": "A" * 20,
    }
    out, requests, _ = fetch_rcsb_fixture(context=context)
    values = [
        json.loads(r.url.params["json"])["query"]["parameters"]["value"]
        for r in requests
        if r.url.host == "search.rcsb.org"
    ]
    assert values == ["P00519", "P00520"]
    assert len([r for r in requests if "/polymer_entity/" in r.url.path]) == 1
    assert out.result.source_context["preferred_accession"] == "P00520"
    assert (
        out.result.records[1].binding_sites[0]["residues"][0]["uniprot_accession"]
        == "P00519"
    )
    assert (
        out.result.records[1].binding_sites[0]["residues"][0]["sequence_check"]
        == "not_checked"
    )


def test_multiple_copies_share_chemistry_and_preserve_sites(
    fetch_rcsb_fixture, rcsb_data
):
    obj = copy.deepcopy(rcsb_data["nonpolymer_entity_instance/1ABC/C"])
    obj["rcsb_id"] = "1ABC.G"
    obj["rcsb_nonpolymer_entity_instance_container_identifiers"]["asym_id"] = "G"
    rcsb_data["nonpolymer_entity_instance/1ABC/G"] = obj
    rcsb_data["nonpolymer_entity/1ABC/2"][
        "rcsb_nonpolymer_entity_container_identifiers"
    ]["asym_ids"].append("G")
    # Duplicate source rows and alternate locations must not multiply residues/sites.
    obj["rcsb_target_neighbors"].append(
        obj["rcsb_target_neighbors"][0] | {"alt_id": "B", "distance": 2.8}
    )
    out, requests, _ = fetch_rcsb_fixture()
    assert out.result.counts["unique_ligands"] == 2
    assert out.result.counts["included_ligand_instances"] == 3
    assert out.result.counts["contact_sites"] == 5
    assert sum("/chemcomp/STI" in r.url.path for r in requests) == 1
    site = next(
        s
        for s in out.result.records[1].binding_sites
        if s["ligand_asym_id"] == "G" and s["target_asym_id"] == "A"
    )
    assert len(site["residues"]) == 1 and site["residues"][0]["distance"] == 2.8


@pytest.mark.parametrize(
    "corruption",
    [
        "count",
        "missing_site",
        "distance",
        "position",
        "chain",
        "model",
        "chemistry",
        "evidence",
        "audit",
        "duplicate",
        "provenance",
        "missing_context",
        "null_complete",
    ],
)
def test_saved_structural_results_reject_corruption(fetch_rcsb_fixture, corruption):
    out, _, _ = fetch_rcsb_fixture()
    data = asdict(out.result)
    record = data["records"][1]
    site = record["binding_sites"][0]
    if corruption == "count":
        data["counts"]["contact_sites"] += 1
    if corruption == "missing_site":
        site.pop("residues")
    if corruption == "distance":
        site["residues"][0]["distance"] = -1
    if corruption == "position":
        site["residues"][0]["uniprot_position"] = True
    if corruption == "chain":
        site["target_asym_id"] = "Z"
    if corruption == "model":
        site["target_model_id"] = True
    if corruption == "chemistry":
        record["structure_provenance"]["ccd_id"] = "CL"
    if corruption == "evidence":
        record["measurement_value"] = "5"
    if corruption == "audit":
        data["selection_audit"]["instances"][0]["included"] = False
    if corruption == "duplicate":
        data["records"].append(copy.deepcopy(record))
    if corruption == "provenance":
        record["structure_provenance"]["status"] = "unavailable"
    if corruption == "missing_context":
        data["source_context"]["entries"] = {}
    if corruption == "null_complete":
        data["counts"]["complete"] = None
    with pytest.raises(ValueError):
        source_result_from_dict(data)


@pytest.mark.parametrize("status", [429, 503])
def test_core_retry_and_throttling(fetch_rcsb_fixture, status):
    calls = []

    def handler(request):
        calls.append(request)
        return (
            httpx.Response(status, headers={"Retry-After": "3"})
            if len(calls) == 1
            else httpx.Response(204)
        )

    out, _, sleeps = fetch_rcsb_fixture(handler=handler)
    assert out.result.status == "no_matches" and len(out.exchanges) == 2
    assert 3 in sleeps and 0.5 in sleeps
    assert out.exchanges[0].request_id == out.exchanges[1].request_id


def test_cli_three_sources_and_rcsb_failure_isolation(
    monkeypatch, tmp_path, payload, empty_uniprot
):
    from typer.testing import CliRunner

    from ligandscout import cli

    original = httpx.Client
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(cli.time, "sleep", lambda _: None)

    def handler(request):
        if request.url.host == "www.bindingdb.org":
            return httpx.Response(200, json=payload)
        if request.url.host == "rest.uniprot.org":
            return httpx.Response(200, json=empty_uniprot)
        assert request.url.host == "search.rcsb.org"
        return httpx.Response(404)

    monkeypatch.setattr(
        cli.httpx,
        "Client",
        lambda **_: original(transport=httpx.MockTransport(handler)),
    )
    result = CliRunner().invoke(cli.app, ["P00699"])
    assert result.exit_code == 1
    assert "Source: BindingDB | Status: success" in result.stdout
    assert "Source: UniProt | Status: no_matches" in result.stdout
    assert "Source: RCSB PDB | Status: failed" in result.stdout
    folder = tmp_path / "results/P00699"
    assert "Contact sites: Unavailable" in result.stdout
    assert "rcsb-ligands.html" not in result.stdout
    assert (folder / "rcsb-ligands.html").exists()
    assert json.loads((folder / "sources/rcsb.json").read_text())["records"] == []


def test_cli_complete_structural_search(
    monkeypatch, tmp_path, payload, empty_uniprot, rcsb_data
):
    from typer.testing import CliRunner

    from ligandscout import cli

    original = httpx.Client
    empty_uniprot["primaryAccession"] = "P00519"
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(cli.time, "sleep", lambda _: None)

    def handler(request):
        if request.url.host == "www.bindingdb.org":
            return httpx.Response(200, json=payload)
        if request.url.host == "rest.uniprot.org":
            return httpx.Response(200, json=empty_uniprot)
        if request.url.host == "search.rcsb.org":
            return httpx.Response(
                200,
                json={
                    "result_type": "polymer_entity",
                    "total_count": 1,
                    "result_set": [{"identifier": "1ABC_1"}],
                },
            )
        assert request.url.host == "data.rcsb.org"
        return httpx.Response(
            200, json=rcsb_data[request.url.path.removeprefix("/rest/v1/core/")]
        )

    monkeypatch.setattr(
        cli.httpx,
        "Client",
        lambda **_: original(transport=httpx.MockTransport(handler)),
    )
    result = CliRunner().invoke(cli.app, ["P00519"])
    assert result.exit_code == 0, result.output
    assert "Source: RCSB PDB | Status: success" in result.stdout
    assert "Contact sites: 3" in result.stdout
    folder = tmp_path / "results/P00519"
    for source in ("bindingdb", "uniprot", "rcsb"):
        data = json.loads((folder / f"sources/{source}.json").read_text())
        assert data["schema_version"] == 5
        assert source_result_from_dict(data).source
    assert "Structural ligand records: 2" in result.stdout


def test_ccd_ids_link_to_ligand_pages(fetch_rcsb_fixture, tmp_path):
    from html.parser import HTMLParser

    class Links(HTMLParser):
        def __init__(self):
            super().__init__()
            self.links = {}
            self.href = None

        def handle_starttag(self, tag, attrs):
            if tag == "a":
                self.href = dict(attrs).get("href")

        def handle_data(self, data):
            if self.href:
                self.links[data] = self.href

        def handle_endtag(self, tag):
            if tag == "a":
                self.href = None

    outcome, _, _ = fetch_rcsb_fixture()
    folder = save_results(outcome, root=tmp_path)
    parser = Links()
    parser.feed((folder / "rcsb-ligands.html").read_text())
    assert outcome.result.records
    for record in outcome.result.records:
        assert parser.links[record.compound_id] == (
            "https://www.rcsb.org/ligand/" + record.compound_id.removeprefix("CCD:")
        )


def test_mismatches_highlighted_without_summary_warnings(fetch_rcsb_fixture, tmp_path):
    from ligandscout.report import render_summary, render_index

    out, _, _ = fetch_rcsb_fixture(
        context={
            "primary_accession": "P00519",
            "secondary_accessions": [],
            "sequence": "GGGGGGGGGG",
        }
    )
    mismatches = [
        w for w in out.result.warnings if w["message"].endswith(": mapped, mismatch.")
    ]
    assert mismatches
    out.result.warnings.append(
        {"code": "other", "message": "Keep this warning", "role": "rcsb_entry"}
    )
    folder = save_results(out, root=tmp_path)
    page = (folder / "rcsb-ligands.html").read_text()
    assert '<li class="residue-mismatch">' in page
    assert "sequence mismatch" in page
    for text in (
        page,
        render_summary(out, folder=folder),
        render_index(out),
        (folder / "report.md").read_text(),
    ):
        assert "mapped, mismatch." not in text
        assert "Keep this warning" in text
    # Detailed evidence and raw warnings remain available in the saved data.
    assert "mapped, mismatch." in (folder / "sources/rcsb.json").read_text()
