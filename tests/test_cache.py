"""Completed searches are reused independently, with their provenance preserved."""

import json

import httpx
import pytest
from typer.testing import CliRunner

from ligandscout.cache import load_completed
from ligandscout.search import run_search
from ligandscout.storage import save_results


@pytest.fixture
def saved(tmp_path, payload, empty_uniprot, now):
    requests = []

    def handler(request):
        requests.append(request.url.host)
        if request.url.host == "search.rcsb.org":
            return httpx.Response(204)
        if request.url.host == "rest.uniprot.org":
            return httpx.Response(200, json=empty_uniprot)
        return httpx.Response(200, json=payload)

    with httpx.Client(transport=httpx.MockTransport(handler)) as client:

        def search(**kwargs):
            return run_search(
                "P00699",
                client=client,
                sleep=lambda _: None,
                now=now,
                root=tmp_path,
                **kwargs,
            )

        outcome = search()
        save_results(outcome, root=tmp_path)
        requests.clear()
        yield tmp_path, search, requests, handler


def test_reuse_preserves_source_files_and_rebuilds_reports(saved, caplog):
    root, search, requests, _ = saved
    folder = root / "P00699"
    paths = list((folder / "sources").rglob("*"))
    snapshot = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in paths if p.is_file()}
    for name in ("index.html", "report.md", "ligands.tsv"):
        (folder / name).unlink()
    for page in folder.glob("*-ligands.html"):
        page.write_text("Old report from a previous renderer")
    with caplog.at_level("INFO"):
        outcome = search()
    assert requests == []
    assert all(s.reused for s in outcome.sources)
    assert "Skipping BindingDB" in caplog.text
    assert [s.result.status for s in outcome.sources] == [
        "success",
        "no_matches",
        "no_matches",
    ]
    save_results(outcome, root=root)
    for page in folder.glob("*-ligands.html"):
        assert "Old report from a previous renderer" not in page.read_text()
        assert "<!doctype html>" in page.read_text().lower()
    assert all(
        (p.read_bytes(), p.stat().st_mtime_ns) == old for p, old in snapshot.items()
    )
    assert all(
        (folder / name).is_file() for name in ("index.html", "report.md", "ligands.tsv")
    )


def test_force_queries_every_source(saved):
    _, search, requests, _ = saved
    outcome = search(force=True)
    assert len(requests) == 3
    assert not any(s.reused for s in outcome.sources)


@pytest.mark.parametrize(
    "relative,host",
    [
        ("sources/bindingdb.json", "www.bindingdb.org"),
        ("bindingdb-ligands.html", "www.bindingdb.org"),
        ("sources/bindingdb.raw.json", "www.bindingdb.org"),
        ("sources/uniprot-requests.json", "rest.uniprot.org"),
        ("sources/uniprot-raw/000001.json", "rest.uniprot.org"),
        ("sources/rcsb-requests.json", "search.rcsb.org"),
    ],
)
def test_missing_artifact_refreshes_only_affected_source(saved, relative, host):
    root, search, requests, _ = saved
    (root / "P00699" / relative).unlink()
    outcome = search()
    assert requests == [host]
    assert sum(s.reused for s in outcome.sources) == 2
    save_results(outcome, root=root)
    assert all(search_source(root, s.result.source) for s in outcome.sources)


def search_source(root, source):
    return load_completed(root, "P00699", source)


@pytest.mark.parametrize(
    "mutation", ["checksum", "path", "accession", "malformed", "failed"]
)
def test_invalid_saved_evidence_is_refetched(saved, mutation):
    root, search, requests, _ = saved
    folder = root / "P00699/sources"
    if mutation in ("checksum", "path"):
        path = folder / "uniprot-requests.json"
        data = json.loads(path.read_text())
        data["exchanges"][0]["sha256" if mutation == "checksum" else "raw_path"] = (
            "invalid"
        )
    else:
        path = folder / "uniprot.json"
        data = json.loads(path.read_text())
        if mutation == "accession":
            data["accession"] = "P00519"
        elif mutation == "failed":
            data["status"] = "failed"
    path.write_text("{" if mutation == "malformed" else json.dumps(data))
    outcome = search()
    assert requests == ["rest.uniprot.org"]
    assert [s.reused for s in outcome.sources] == [True, False, True]


def test_cli_force_and_default(saved, monkeypatch):
    from ligandscout import cli

    root, _, requests, handler = saved
    # CLI uses cwd/results; move the fixture output under that location.
    results = root / "results"
    results.mkdir()
    (root / "P00699").rename(results / "P00699")
    monkeypatch.chdir(root)
    real_client = httpx.Client
    monkeypatch.setattr(
        cli.httpx,
        "Client",
        lambda **kw: real_client(transport=httpx.MockTransport(handler)),
    )
    runner = CliRunner()
    reused = runner.invoke(cli.app, ["P00699"])
    assert reused.exit_code == 0, reused.output
    assert requests == []
    assert "Reused saved results" in reused.stdout
    forced = runner.invoke(cli.app, ["P00699", "--force"])
    assert forced.exit_code == 0, forced.output
    assert len(requests) == 3
    assert "Reused saved results" not in forced.stdout
