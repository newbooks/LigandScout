"""Reports, saved provenance, and overwrite failures."""

import csv
import io
import json
from dataclasses import asdict, replace

import httpx
import pytest

from ligandscout.models import SearchOutcome, SourceError, upgrade_result
from ligandscout.report import COLUMNS, render_report, render_tsv
from ligandscout.storage import save_results


def test_report_and_schema(positive, tmp_path):
    folder = save_results(positive, root=tmp_path)
    result = upgrade_result(positive.result)
    data = json.loads((folder / "sources/bindingdb.json").read_text())
    assert data == asdict(result)
    assert set(data) == {
        "schema_version",
        "accession",
        "source",
        "retrieved_at",
        "request",
        "source_version",
        "status",
        "records",
        "error",
        "policy",
        "warnings",
        "counts",
        "selection_audit",
        "source_context",
    }
    assert (folder / "sources/bindingdb.raw.json").read_bytes() == positive.raw_body
    assert (folder / "report.md").read_text() == render_report(result)
    report = render_report(result)
    assert "Unique compounds: 1" in report and "Measurements: 1" in report
    assert "unspecified" in report and "pubmed.ncbi.nlm.nih.gov/16248836/" in report
    assert "doi.org/10.2174/" in report
    assert (folder / "ligands.tsv").read_text() == render_tsv(result)
    rows = list(csv.DictReader(io.StringIO(render_tsv(result)), delimiter="\t"))
    assert rows[0]["measurement_unit"] == ""
    assert rows[0]["source"] == "BindingDB"
    assert rows[0]["retrieved_at"] == result.retrieved_at


def test_tsv_quotes_and_markdown_escapes(positive):
    row = positive.result.records[0]
    row.smiles = "C\tC\nC"
    row.measurement_value = "> 1|[x]*\n<script>"
    row.pmid = "https://evil.invalid/"
    row.doi = "10.123/a)b|x"
    tsv = render_tsv(positive.result)
    parsed = list(csv.DictReader(io.StringIO(tsv, newline=""), delimiter="\t"))
    assert parsed[0]["smiles"] == row.smiles
    assert parsed[0]["measurement_value"] == row.measurement_value
    report = render_report(positive.result)
    assert "&gt; 1\\|\\[x\\]\\* &lt;script&gt;" in report
    assert "%29b%7Cx" in report
    assert "(https://evil.invalid/)" not in report


def test_duplicate_counts_and_unknown_reference(positive):
    positive.result.records *= 2
    positive.result.records[0].doi = "unknown"
    report = render_report(positive.result)
    assert "Unique compounds: 1" in report and "Measurements: 2" in report
    assert "unknown" in report and "doi.org" not in report


@pytest.mark.parametrize(
    "body,extension",
    [(b"", "txt"), (b"<html>bad</html>", "txt"), (b'""', "json"), (b"NaN", "txt")],
)
def test_raw_variant_and_no_stale_records(positive, fetch, tmp_path, body, extension):
    folder = save_results(positive, root=tmp_path)
    (folder / "keep.txt").write_text("user data")
    (folder / "sources/other.json").write_text("other source")
    next_outcome, _ = fetch(lambda request: httpx.Response(200, content=body))
    save_results(next_outcome, root=tmp_path)
    assert (folder / f"sources/bindingdb.raw.{extension}").read_bytes() == body
    other = "json" if extension == "txt" else "txt"
    assert not (folder / f"sources/bindingdb.raw.{other}").exists()
    rows = list(
        csv.reader(io.StringIO((folder / "ligands.tsv").read_text()), delimiter="\t")
    )
    assert rows == [COLUMNS]
    assert json.loads((folder / "sources/bindingdb.json").read_text())["records"] == []
    assert (folder / "keep.txt").read_text() == "user data"
    assert (folder / "sources/other.json").read_text() == "other source"
    assert not list(folder.rglob(".ligandscout-*"))


def test_transport_failure_removes_both_raw_variants(positive, tmp_path):
    folder = save_results(positive, root=tmp_path)
    (folder / "sources/bindingdb.raw.txt").write_text("stale")
    failed = replace(
        positive.result,
        status="failed",
        records=[],
        error=SourceError("transport", "Unavailable", None),
    )
    save_results(SearchOutcome(failed, None), root=tmp_path)
    assert not list((folder / "sources").glob("*.raw.*"))
    report = (folder / "report.md").read_text()
    assert "Search failed" in report and "No BindingDB records" not in report


def test_publish_failure_cleans_temps(positive, tmp_path, monkeypatch):
    from ligandscout import storage

    real_replace = storage.os.replace

    def fail_report(source, destination):
        if destination.name == "report.md":
            raise OSError("disk full")
        real_replace(source, destination)

    monkeypatch.setattr(storage.os, "replace", fail_report)
    with pytest.raises(OSError):
        save_results(positive, root=tmp_path)
    assert not list(tmp_path.rglob(".ligandscout-*"))
    assert not (tmp_path / "P00699/sources/bindingdb.json").exists()


def test_staging_failure_preserves_old_files(positive, tmp_path, monkeypatch):
    from ligandscout import storage

    folder = save_results(positive, root=tmp_path)
    before = {p: p.read_bytes() for p in folder.rglob("*") if p.is_file()}
    monkeypatch.setattr(
        storage.tempfile,
        "NamedTemporaryFile",
        lambda **kwargs: (_ for _ in ()).throw(OSError("disk full")),
    )
    with pytest.raises(OSError):
        save_results(positive, root=tmp_path)
    assert all(path.read_bytes() == body for path, body in before.items())


def test_interrupt_cleans_temps(positive, tmp_path, monkeypatch):
    from ligandscout import storage

    monkeypatch.setattr(
        storage.os, "replace", lambda *args: (_ for _ in ()).throw(KeyboardInterrupt())
    )
    with pytest.raises(KeyboardInterrupt):
        save_results(positive, root=tmp_path)
    assert not list(tmp_path.rglob(".ligandscout-*"))


def test_summary_is_bounded_without_dropping_saved_records(positive, tmp_path):
    from ligandscout.report import render_summary

    original = positive.result.records[0]
    positive.result.records = [
        replace(
            original,
            compound_id=str(index),
            target_description=f"Target {index}",
            measurement_type=f"Type {index}",
        )
        for index in range(1000)
    ]
    summary = render_summary(positive.result, folder=tmp_path)
    assert len(summary.splitlines()) <= 23
    assert (
        f"Unique ligands listed in BindingDB for {positive.result.accession}: 1000"
        in summary
    )
    assert "Measurements: 1000" in summary
    assert summary.count("995 additional labels") == 2
    assert "not confirmed binding partners" in summary
    assert "have not been independently verified" in summary
    assert (
        len(
            list(
                csv.DictReader(io.StringIO(render_tsv(positive.result)), delimiter="\t")
            )
        )
        == 1000
    )
    report = render_report(positive.result)
    assert report.count("?monomerid=") == 1000
    assert "Reported target" in report and "Target 999" in report


def test_summary_counts_labels_and_sanitizes_them(positive, tmp_path):
    from ligandscout.report import render_summary

    row = positive.result.records[0]
    positive.result.records = [replace(row, target_description="Major")] * 3 + [
        replace(
            row,
            target_description="Other\n\x1btarget",
            measurement_type="Other\n\x1btype",
        )
    ]
    summary = render_summary(positive.result, folder=tmp_path)
    assert "Reported targets (measurement counts):" in summary
    assert "Major: 3" in summary and "Other  target: 1" in summary
    assert summary.index("Major: 3") < summary.index("Other  target: 1")
    assert (
        f"Unique ligands listed in BindingDB for {positive.result.accession}: 1"
        in summary
    )
    assert "Measurements: 4" in summary
    assert "Ki: 3" in summary and "Other  type: 1" in summary
    assert "\x1b" not in summary
    assert summary.index("Ki: 3") < summary.index("Other  type: 1")
    assert summary.endswith(
        "Full measurements: ligands.tsv | Details: index.html | Markdown: report.md\n"
    )
    assert "bindingdb-ligands.html" not in summary


def test_html_unique_ligands_and_summary_link(positive, tmp_path):
    from ligandscout.report import render_html

    row = positive.result.records[0]
    positive.result.records = [row, replace(row, target_description="Other target")]
    folder = save_results(positive, root=tmp_path)
    page = (folder / "bindingdb-ligands.html").read_text()
    assert page == render_html(positive.result)
    assert "Unique ligands: 1 | Measurements: 2" in page
    assert page.count("<tr><td><a ") == 1
    assert "2 measurements" in page and "Other target" in page
    report = (folder / "report.md").read_text()
    assert "| Database | Ligands Count | Ligands Report |" in report
    assert "| BindingDB | 1 | [View ligands](bindingdb-ligands.html) |" in report


def test_html_escapes_remote_text_and_builds_safe_links(positive):
    from ligandscout.report import render_html

    row = positive.result.records[0]
    row.compound_name = '<script>alert("x")</script>'
    row.target_description = "<img src=x onerror=alert(1)>"
    row.measurement_value = "<bad>"
    row.source_url = "javascript:alert(1)"
    page = render_html(positive.result)
    assert "<script>" not in page and "<img" not in page
    assert "&lt;script&gt;" in page and "&lt;bad&gt;" in page
    assert "javascript:" not in page
    assert "https://www.bindingdb.org/" in page


@pytest.mark.parametrize("status", ["failed", "no_matches"])
def test_html_empty_results_replace_previous_report(positive, tmp_path, status):
    folder = save_results(positive, root=tmp_path)
    result = replace(
        positive.result,
        records=[],
        status=status,
        error=SourceError("transport", "Service unavailable", None)
        if status == "failed"
        else None,
    )
    save_results(SearchOutcome(result, None), root=tmp_path)
    page = (folder / "bindingdb-ligands.html").read_text()
    count = "Unavailable" if status == "failed" else "0"
    assert f"Unique ligands: {count} | Measurements: {count}" in page
    assert "<tr><td><a " not in page
    assert ("Search failed:" in page) == (status == "failed")
    assert ("No records found." in page) == (status == "no_matches")
    if status == "failed":
        assert "| BindingDB | Unavailable |" in (folder / "report.md").read_text()


def test_html_smiles_visible_unique_and_escaped(positive):
    import re
    from html import unescape

    from ligandscout.report import render_html

    row = positive.result.records[0]
    first = 'C/C=C\\C.[Na+]&<tag>"'
    second = "C" * 500
    result = replace(
        positive.result,
        source="BindingDB",
        records=[
            replace(row, smiles=first),
            replace(row, smiles=first),
            replace(row, smiles=second),
            replace(row, smiles=None),
        ],
    )
    page = render_html(result)
    assert "<th>SMILES</th>" in page
    blocks = re.findall(r'<code class="smiles">(.*?)</code>', page, re.DOTALL)
    assert [unescape(value) for value in blocks] == sorted([first, second])
    assert "<tag>" not in page
    assert page.index('<code class="smiles">') < page.index("<details>")
    assert "user-select:all" in page
    assert "overflow-wrap:anywhere" in page


def test_html_missing_smiles_is_explicit(positive):
    from ligandscout.report import render_html

    row = positive.result.records[0]
    result = replace(
        positive.result,
        records=[
            replace(row, smiles=None),
            replace(row, smiles="   "),
        ],
    )
    page = render_html(result)
    assert "<th>SMILES</th>" in page
    assert '<code class="smiles">' not in page
    assert "<td>Not supplied</td>" in page


def test_old_source_artifacts_are_not_loaded_or_linked(positive, tmp_path):
    folder = tmp_path / positive.result.accession
    sources = folder / "sources"
    sources.mkdir(parents=True)
    old_files = {
        sources / "chembl.json": b'{"records": "old"}',
        sources / "pubchem.json": b'{"records": "old"}',
        folder / "chembl-ligands.html": b"old page",
        folder / "pubchem-ligands.html": b"old page",
    }
    for path, body in old_files.items():
        path.write_bytes(body)
    save_results(positive, root=tmp_path)
    for path, body in old_files.items():
        assert path.read_bytes() == body
    for name in ("report.md", "bindingdb-ligands.html", "ligands.tsv"):
        text = (folder / name).read_text()
        assert "chembl" not in text.lower() and "pubchem" not in text.lower()
    rows = list(
        csv.DictReader(
            io.StringIO((folder / "ligands.tsv").read_text()), delimiter="\t"
        )
    )
    assert {row["source"] for row in rows} == {"BindingDB"}


@pytest.mark.parametrize("version", [1, 2, 3])
def test_bindingdb_saved_schema_compatibility(positive, version):
    from ligandscout.models import source_result_from_dict

    current = upgrade_result(positive.result)
    data = asdict(positive.result if version == 1 else current)
    data["schema_version"] = version
    if version in (2, 3):
        from ligandscout.models import ANNOTATION_COUNTS

        for key in ("record_kind", *ANNOTATION_COUNTS):
            data["counts"].pop(key)
    if version == 2:
        for key in ("candidate_assays", "eligible_assays"):
            data["counts"].pop(key)
        data["selection_audit"].pop("assays")
    assert source_result_from_dict(data) == current


def test_removed_source_results_are_rejected(positive):
    from ligandscout.models import source_result_from_dict

    data = asdict(upgrade_result(positive.result))
    data["source"] = "RemovedSource"
    with pytest.raises(ValueError, match="Unknown source"):
        source_result_from_dict(data)


@pytest.mark.parametrize("status", ["success", "no_matches", "failed"])
def test_index_summary_links_and_rerun_overwrite(positive, tmp_path, status):
    from html.parser import HTMLParser

    from ligandscout.report import render_index

    class Links(HTMLParser):
        def __init__(self):
            super().__init__()
            self.hrefs = []

        def handle_starttag(self, tag, attrs):
            if tag == "a":
                self.hrefs.append(dict(attrs)["href"])

    folder = save_results(positive, root=tmp_path)
    result = replace(
        positive.result,
        status=status,
        records=positive.result.records if status == "success" else [],
        error=SourceError("http", "Service unavailable", 503)
        if status == "failed"
        else None,
    )
    save_results(SearchOutcome(result, positive.raw_body), root=tmp_path)
    page = (folder / "index.html").read_text()
    assert page == render_index(result)
    assert f"<td>{status}</td>" in page
    count = "Unavailable" if status == "failed" else "1" if status == "success" else "0"
    assert f"<td>{count}</td><td>{count}</td>" in page
    assert result.retrieved_at in page
    assert "not confirmed binding partners" in page
    assert ("search failed:" in page) == (status == "failed")
    assert ("No BindingDB records found" in page) == (status == "no_matches")
    links = Links()
    links.feed(page)
    assert set(links.hrefs) == {
        "bindingdb-ligands.html",
        "sources/bindingdb.json",
        "ligands.tsv",
        "report.md",
    }
    assert all((folder / href).is_file() for href in links.hrefs)
    assert 'href="index.html"' in (folder / "bindingdb-ligands.html").read_text()
    assert "chembl" not in page.lower() and "pubchem" not in page.lower()


def test_index_escapes_errors_and_warnings(positive):
    from ligandscout.report import render_index

    result = replace(
        positive.result,
        records=[],
        status="failed",
        error=SourceError("schema", "<script>alert('bad')</script>", 200),
        warnings=[
            {"code": "test", "message": '<img src=x onerror="bad">', "role": None}
        ],
    )
    page = render_index(result)
    assert "<script>" not in page and "<img" not in page
    assert "&lt;script&gt;" in page and "&lt;img" in page
