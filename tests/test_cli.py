"""Check installed commands and mocked end-to-end searches."""

import subprocess
import sys

import httpx
import pytest
from typer.testing import CliRunner

from ligandscout.cli import app

runner = CliRunner()


def test_help_aliases():
    short = runner.invoke(app, ["-h"])
    long = runner.invoke(app, ["--help"])
    assert short.exit_code == long.exit_code == 0
    assert short.stdout == long.stdout
    assert "ACCESSION_ID" in long.stdout
    assert "--simple" not in long.stdout
    assert "--explore" not in long.stdout


def test_search_queries_and_saves(monkeypatch, tmp_path, payload, empty_uniprot):
    from ligandscout import cli

    monkeypatch.chdir(tmp_path)
    real_client = httpx.Client

    def handler(request):
        if request.url.host == "search.rcsb.org":
            return httpx.Response(204)
        if request.url.host == "rest.uniprot.org":
            return httpx.Response(200, json=empty_uniprot)
        assert request.url.host == "www.bindingdb.org"
        assert request.url.params["uniprot"] == "P00699"
        return httpx.Response(200, json=payload)

    monkeypatch.setattr(
        cli.httpx,
        "Client",
        lambda **kwargs: real_client(transport=httpx.MockTransport(handler)),
    )
    result = runner.invoke(app, [" p00699 "])
    assert result.exit_code == 0, result.output
    folder = tmp_path / "results/P00699"
    assert "Full measurements: ligands.tsv" in result.stdout
    assert "| Compound |" not in result.stdout
    assert "| Compound |" in (folder / "report.md").read_text()
    assert "Status: success" in result.stdout
    assert "[INFO]: Searching BindingDB for P00699" in result.stderr
    assert "[INFO]" not in result.stdout


@pytest.mark.parametrize("flag", ["--simple", "--explore"])
@pytest.mark.parametrize("accession", [[], ["P00699"]])
def test_removed_flags_are_rejected_without_side_effects(
    flag, accession, monkeypatch, tmp_path
):
    from ligandscout import cli

    monkeypatch.chdir(tmp_path)

    def forbidden(**kwargs):
        pytest.fail("Removed flags must not create a network client")

    monkeypatch.setattr(cli.httpx, "Client", forbidden)
    result = runner.invoke(app, [*accession, flag])
    assert result.exit_code == 2
    assert "No such option" in result.output
    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize("args", [[], ["--log-level", "DEBUG"]])
def test_missing_accession_displays_help(args):
    result = runner.invoke(app, args)
    assert result.exit_code == 0
    assert result.stdout == runner.invoke(app, ["--help"]).stdout
    assert result.stderr == ""


@pytest.mark.parametrize("args", [["P00699", "--unknown"], ["P00699", "-l", "INVALID"]])
def test_invalid_arguments(args):
    assert runner.invoke(app, args).exit_code == 2


@pytest.mark.parametrize(
    "command", [["ligandscout"], [sys.executable, "-m", "ligandscout"]]
)
def test_installed_entry_points(command):
    result = subprocess.run(
        [*command, "--help"], capture_output=True, text=True, check=False
    )
    assert result.returncode == 0, result.stderr
    assert "ACCESSION_ID" in result.stdout

    no_args = subprocess.run(command, capture_output=True, text=True, check=False)
    assert no_args.returncode == 0, no_args.stderr
    assert no_args.stdout == result.stdout
    assert no_args.stderr == ""


@pytest.mark.parametrize(
    "args",
    [
        ["LYZ"],
        ["P00699-2"],
        ["../P00699"],
        ["--explore", "P00699"],
        ["--help"],
        [],
        ["--unknown"],
        ["--log-level"],
    ],
)
def test_no_network_or_files_for_nonsearch_commands(args, monkeypatch, tmp_path):
    from ligandscout import cli

    monkeypatch.chdir(tmp_path)

    def forbidden(**kwargs):
        pytest.fail("Network client should not be created")

    monkeypatch.setattr(cli.httpx, "Client", forbidden)
    runner.invoke(app, args)
    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize(
    "status,body,expected",
    [
        (200, b'{"getLindsByUniprotsResponse":{"affinities":[]}}', 0),
        (200, b"<html>bad</html>", 1),
        (404, b"", 1),
        (200, b"", 1),
    ],
)
def test_cli_empty_and_failed_results(
    status, body, expected, monkeypatch, tmp_path, empty_uniprot
):
    from ligandscout import cli

    monkeypatch.setattr(cli.time, "sleep", lambda seconds: None)

    monkeypatch.chdir(tmp_path)
    real_client = httpx.Client

    def handler(request):
        if request.url.host == "search.rcsb.org":
            return httpx.Response(204)
        if request.url.host == "rest.uniprot.org":
            return httpx.Response(200, json=empty_uniprot)
        assert request.url.host == "www.bindingdb.org"
        return httpx.Response(status, content=body)

    monkeypatch.setattr(
        cli.httpx,
        "Client",
        lambda **kwargs: real_client(transport=httpx.MockTransport(handler)),
    )
    result = runner.invoke(app, ["P00699"])
    assert result.exit_code == expected, result.output
    assert "Full measurements: ligands.tsv" in result.stdout
    assert (tmp_path / "results/P00699/report.md").exists()
    if expected:
        assert (
            "Search failed" in result.stdout
            and "No BindingDB records" not in result.stdout
        )
        assert "[ERROR]:" in result.stderr
    else:
        assert "No BindingDB records found for P00699" in result.stdout


def test_cli_output_failure(monkeypatch, tmp_path, positive):
    from ligandscout import cli

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(cli, "run_search", lambda *args, **kwargs: positive)
    monkeypatch.setattr(
        cli,
        "save_results",
        lambda *args, **kwargs: (_ for _ in ()).throw(OSError("disk full")),
    )
    result = runner.invoke(app, ["P00699"])
    assert result.exit_code == 1 and result.stdout == ""
    assert "incomplete or mixed" in result.stderr


def test_cli_interrupt(monkeypatch, tmp_path):
    from ligandscout import cli

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        cli,
        "run_search",
        lambda *args, **kwargs: (_ for _ in ()).throw(KeyboardInterrupt()),
    )
    result = runner.invoke(app, ["P00699"])
    assert result.exit_code == 130
    assert "Interrupted" in result.stderr
