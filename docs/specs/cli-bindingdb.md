# CLI and BindingDB implementation specification

Status: implemented. The selected sources are BindingDB, UniProtKB/Swiss-Prot, and RCSB PDB. The command queries BindingDB and UniProt; RCSB CCD supplies optional exact-ID chemical enrichment. Experimental PDB interaction search is not implemented.

`pip install -e .` installs the Typer command and `python -m ligandscout` entry point. The [UniProt integration specification](uniprot-integration.md) defines current two-source orchestration, schema v4, request auditing, annotation/site reporting, and optional CCD enrichment. It supersedes the original single-source orchestration and schema descriptions below. BindingDB's request/parser/retry/measurement behavior is unchanged.

## Original BindingDB-only orchestration and saved-format contract (historical)

- `run_search(accession, *, client, sleep, now) -> RunOutcome` normalizes the accession, fetches BindingDB once (subject to adapter retries), upgrades its metadata, and returns one source. It performs no other service requests.
- `RunOutcome.sources` remains a list to support future selected adapters. Rendering and storage retain this interface, but currently accept BindingDB source results only.
- BindingDB's adapter produces its internal version-1 result. Orchestration and storage upgrade it to the existing shared schema version 3, retaining optional measurement/assay fields, `policy`, `warnings`, `counts`, and `selection_audit`. BindingDB v1/v2/v3 saved-result reading remains supported; removed source types are rejected. Keeping the schema avoids needless changes to existing BindingDB outputs.
- Save `sources/bindingdb.json`, the final raw response, `bindingdb-ligands.html`, `ligands.tsv`, `report.md`, and `index.html`. HTML provides one row per unique compound with selectable SMILES and expandable measurement details. Markdown contains a database/count/HTML-link table plus full measurements. Stdout remains compact.
- Stage all artifacts before publishing. Replace raw files, database HTML, TSV, Markdown, index HTML, then normalized JSON; remove obsolete BindingDB raw variants after publication. Atomicity remains per file.
- Regenerate reports only from current-run BindingDB records. Preserve historical files from removed integrations and unrelated user files without reading or linking them. Do not recursively delete existing result directories.
- A BindingDB failure is saved and produces exit code 1; a completed empty search exits 0. UniProt and PDB are not included as queried or failed sources until their adapters exist.

The sections below preserve the original BindingDB data contract (including its internal v1 fields). The current shared schema/output extensions above take precedence over older single-source signatures or v1 serialization descriptions.

## HTML summary index

Current two-source indexes use “Records” for the shared count column: BindingDB activity measurements versus UniProt ligand-annotation groups. Source-specific notes identify the evidence differences.

Every saved search creates or overwrites `results/ACCESSION/index.html` alongside the existing reports. Render it from the current `RunOutcome` through `render_index(value) -> str`; it must not load prior source files or initiate network requests. The page works locally without external scripts or stylesheets.

Show the accession and the evidence disclaimer. Include one table row per queried database with its status, unique ligand count, measurement count, UTC retrieval time, relative link to `SLUG-ligands.html`, and relative link to `sources/SLUG.json`. Failed sources show `Unavailable` counts and escaped failure messages; completed empty searches show zero and an explicit no-match explanation. Include escaped source warnings and links to `ligands.tsv` and `report.md`. Do not list planned adapters as queried sources or add a cross-database total.

Database result pages link back to `index.html`. The terminal summary identifies the index as the entry page and lists only `ligands.tsv`, `index.html`, and `report.md`; individual database HTML filenames are linked from the index. Stage and overwrite the index through the existing output publication path, including no-match and failed-source runs; output errors retain the existing exit behavior. Test generated links against files saved in a temporary accession directory, success/empty/failure counts, rerun overwrites, and HTML escaping of source errors/warnings.

## BindingDB implementation workflow

1. Parse the command and handle help and invalid options without network or filesystem changes.
2. Normalize and validate the supplied accession. Use the same normalized value for the API request, saved metadata, and result directory.
3. Request BindingDB data for that accession through the API, applying the retry policy below. Do not require a preliminary UniProt call for this integration.
4. Preserve the final HTTP response body, validate the response, and map each measurement to the shared record schema. A transport failure has no response body; do not reuse a body from an earlier attempt.
5. Build a source result with status `success`, `no_matches`, or `failed`. Retain all measurements; group by compound ID only when calculating the unique-ligand count.
6. Save the source result and available raw response in `results/ACCESSION_ID/sources/`, overwriting the previous BindingDB output.
7. Generate `ligands.tsv` and `report.md` from this invocation's normalized results, then print a compact terminal summary and return the appropriate exit code. A failure report must label the source as failed rather than claiming that no ligands were found.

Future UniProt and RCSB PDB adapters will contribute through the shared result structure. Other database integrations are out of scope.

## Command contract

```text
ligandscout ACCESSION_ID [--log-level LEVEL]
```

```bash
ligandscout P00699
ligandscout P00176 --log-level DEBUG
```

| Argument or option | Behavior |
| --- | --- |
| `ACCESSION_ID` | Positional UniProt accession required to search. Accept one accession per invocation; if omitted, display help and exit successfully. |
| `--log-level`, `-l` | Default `INFO`; accept `DEBUG`, `INFO`, `WARNING`, `ERROR`, and `CRITICAL`, case-insensitively. Reject other values. |
| `--help`, `-h` | Show matching help text and exit successfully without network access. |

If the accession is omitted, display full help on stdout and exit with code 0 before starting a search. The accession invokes the sole search workflow; no mode selector or `search` subcommand is required.

Strip surrounding whitespace and uppercase the accession. Validate canonical six- or ten-character UniProt accession syntax before requesting data. This first version assumes the accession belongs to Swiss-Prot; it does not verify reviewed status or resolve obsolete accessions. Reject names, FASTA sequences, accession lists, and isoform suffixes with a clear usage error. Syntax validation alone does not prove that an entry exists.

Use `re.fullmatch` with this pattern after normalization:

```text
[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}
```

Examples: `p00699` and ` P00699 ` normalize to `P00699`; `A0A024RBG1` is syntactically valid; `P00699-2`, `LYZ`, an empty string, and `../P00699` are invalid. A syntactically valid accession with no BindingDB coverage is `no_matches`, not an accession-validation failure.

Validation order is: Typer parsing and explicit help; missing accession help; log-level validation; accession validation; search. Unknown options and missing option values remain Typer usage errors even without an accession. Keep the positional argument optional in the Typer signature to support missing-accession help.

## Typer and logging conventions

Use a `src/ligandscout/` package and a console entry point `ligandscout = "ligandscout.cli:main"`. Support `python -m ligandscout` through the same entry point.

Keep the Typer layer thin: declare arguments/options, configure logging, invoke the search service, render results, and translate expected exceptions into exit codes. Use `typer.Argument`, `typer.Option`, `typer.echo`, `typer.BadParameter`, and `typer.Exit`. Ensure Typer exposes a single root command rather than interpreting the accession as a subcommand.

Reference conventions inspected in the sibling repositories:

- `../PClass/src/pclass/commands/production.py`: explicit Typer options, `-h`/`--help`, and `--log-level`/`-l`.
- `../PClass/src/pclass/commands/common.py`: standard Python logging with `[%(levelname)s]: %(message)s`.
- `../PClass-Dark/src/pclass_dark/cli.py`: a `main()` wrapper, library operations outside the CLI, and concise exception handling.
- `../PClass-Dark/tests/test_cli.py`: Typer `CliRunner`, entry-point checks, and logging assertions.

These paths are repository references relative to the project root, not runtime dependencies. Adopt PClass's colon punctuation consistently. Use module loggers from `logging.getLogger(__name__)`; logs go to stderr and the compact summary goes to stdout. Log the accession, source, retries, and result counts. Expected service failures produce concise errors without tracebacks at the default log level; debugging information may be logged at `DEBUG`.

## BindingDB request

Use HTTPS and the documented UniProt endpoint:

```text
GET https://www.bindingdb.org/rest/getLigandsByUniprots
    ?uniprot=ACCESSION_ID
    &response=application/json
```

No registration, API key, sequence, or local database download is required. Omit the optional affinity cutoff to avoid intentionally excluding weak measurements. Do not use a sequence-similarity, PDB, protein-name, or homolog fallback when there are no matches.

Use `httpx.Client` and add `httpx>=0.27,<1` to runtime dependencies when implementing. Pass query parameters through `params`, not manual URL concatenation. Send `Accept: application/json` and `User-Agent: ligandscout/0.1.0` (derive the version from installed package metadata). Keep TLS verification enabled. Handle redirects explicitly with `follow_redirects=False`: resolve each Location against the current URL, allow up to five HTTPS redirects within one attempt, and reject missing/invalid Location values or HTTP destinations. Record redirect-policy failures as kind `http`, with the redirect response's status/body when available. No authentication configuration is needed.

The [BindingDB API documentation](https://www.bindingdb.org/rwd/bind/BindingDBRESTfulAPI.jsp) describes accession lookup and an empty-string no-match response. However, live P0DTC1 searches returned an empty HTTP 200 body followed by populated responses for the identical query. Treat empty bodies and JSON empty strings as ambiguous, retry them, and fail if they persist; only a valid structured empty affinity list establishes no matches. Live responses inspected on November 4, 2026 also used:

```json
{
  "getLindsByUniprotsResponse": {
    "affinities": [
      {
        "query": "Cytochrome P450 2B1",
        "monomerid": "50027784",
        "smile": "C#Cc1cc2ccccc2c2ccccc12",
        "affinity_type": "Ki",
        "affinity": "140",
        "pmid": "16248836",
        "doi": "10.2174/138920005774330639"
      }
    ]
  }
}
```

Preserve the actual API spelling `getLindsByUniprotsResponse` and `smile`. The `query` value is a target description, not an echoed accession. Associate records with the single requested accession through request provenance. The endpoint provides database associations; it does not independently establish target specificity or direct binding for every assay.

Apply this response decision table after retries finish:

| Final response | Result |
| --- | --- |
| No complete response body because the final attempt failed in transport | `failed`, error kind `transport`, no raw file |
| HTTP status other than 200 | `failed`, error kind `http`; retain any complete body |
| HTTP 200 with an empty/ASCII-whitespace-only body, or a decoded JSON empty string | Retry within the shared three-attempt budget; if persistent, `failed` with kind `schema`, status 200, and an explicit indeterminate-results message |
| HTTP 200 with the expected envelope and `affinities: []` | `no_matches` |
| HTTP 200 with the expected envelope and a nonempty list of valid records | `success` |
| HTTP 200 with HTML, malformed JSON, wrong envelope/types, or an invalid record | `failed`, error kind `schema`; retain the complete body |

Inspect the body rather than requiring a particular Content-Type, since services can mislabel JSON. Parse strict JSON: reject nonstandard NaN/Infinity constants. Except for the empty-string case, the JSON root and `getLindsByUniprotsResponse` must be objects, and `affinities` must be a list of objects. Ignore unknown extra fields. Required `monomerid`, `affinity_type`, and `affinity` must be nonempty strings after a whitespace check; require `monomerid` to contain only ASCII digits. Preserve their original strings rather than stripping qualifiers or converting numeric values. Numeric JSON values in these fields are schema errors, not silently coerced.

Optional `query`, `smile`, `pmid`, and `doi` may be absent, null, or strings. Normalize empty/whitespace-only optional strings to null and otherwise preserve them. Other types are schema errors. An invalid row fails the entire source result and produces `records: []`; never silently skip rows or return a partial success.

The endpoint documents no pagination mechanism. Consume the full response without a client-side result cap. Do not claim exhaustive database coverage; if the service later documents pagination or truncation, update this specification and the adapter.

## Data and output

### Shared source-result schema

Use `schema_version: 1`. All fields below are present; unavailable values are JSON null, except `records`, which is always a list.

| Field | Type | Meaning |
| --- | --- | --- |
| `schema_version` | integer | Saved-result schema version, independent of the source database release |
| `accession` | string | Normalized query accession |
| `source` | string | `BindingDB` for this adapter |
| `retrieved_at` | string | UTC ISO 8601 timestamp when the final attempt completed, including on failure |
| `request` | object | `url` and `parameters` used for the source query |
| `source_version` | string or null | Database release if explicitly exposed by the API |
| `status` | string | `success`, `no_matches`, or `failed` |
| `records` | list of objects | One normalized record per returned measurement; empty on no matches or failure |
| `error` | object or null | Failure `kind`, concise `message`, and nullable `http_status`; null otherwise |

Failure kinds are `transport`, `http`, or `schema`. Output-writing errors are reported through logs and exit code 1 because the result itself may not be writable. Do not infer a database release from retrieval time.

Each measurement uses these fields:

| Field | BindingDB mapping |
| --- | --- |
| `compound_id` | `monomerid`, retained as a string |
| `compound_name` | null in this adapter version |
| `smiles` | `smile`, preserved verbatim |
| `inchikey` | null in this adapter version |
| `target_description` | `query`, or null if absent |
| `evidence_type` | `reported_activity`; preserve the measurement type without asserting direct binding |
| `measurement_type` | `affinity_type`, preserved verbatim |
| `measurement_value` | `affinity` as original text, including any qualifier |
| `measurement_unit` | null in this adapter version; display `unspecified` in the report |
| `pmid`, `doi` | Supplied publication identifiers, or null |
| `source_url` | `https://www.bindingdb.org/rwd/bind/chemsearch/marvin/MolStructure.jsp?monomerid=COMPOUND_ID` |

Do not parse chemical structures or query another source to fill missing fields in this increment. Preserve duplicate returned rows, since this endpoint does not expose a stable measurement identifier that would establish whether two rows are the same experiment.

All measurement fields are strings or null; `compound_id`, `evidence_type`, `measurement_type`, `measurement_value`, and `source_url` are never null. The endpoint's published description identifies nM for its affinity cutoff but does not explicitly establish the unit of every returned value. Do not infer returned units from the cutoff parameter. A later change to populate units requires a verified source contract and a specification update.

`request` is exactly `{"url": "https://www.bindingdb.org/rest/getLigandsByUniprots", "parameters": {"uniprot": "ACCESSION_ID", "response": "application/json"}}`, with the accession substituted. `source_version` is null for this API version. `error.http_status` is the final HTTP status for HTTP/schema failures, or null for transport failures. Timestamp format is UTC with a `Z` suffix and microseconds, produced by an injectable clock after the final attempt. Error messages must not include the raw body. JSON files use UTF-8, two-space indentation, and a trailing newline.

### Compiled table contract

Write UTF-8 TSV with this fixed column order:

```text
accession, source, compound_id, compound_name, smiles, inchikey,
target_description, evidence_type, measurement_type, measurement_value,
measurement_unit, pmid, doi, source_url, retrieved_at
```

The commas above delimit the column names for readability; the actual file uses tabs. Repeat source-level accession, source, and retrieval time on each row. Retain the existing appended shared metadata columns in the current TSV format. Encode nulls as empty cells, write headers even with zero records, and use Python's `csv.writer` with delimiter `\t`, `QUOTE_MINIMAL`, and line terminator `\n` (open files with `newline=""`).

Keep a measurement record for each returned row, including:

- Requested accession, source (`BindingDB`), and source target description.
- BindingDB monomer ID as a string, supplied SMILES, affinity type, and original affinity text.
- PMID and DOI when supplied, plus a compound source link constructed from the monomer ID.
- Request URL/parameters, UTC retrieval time, and database version if exposed; otherwise version is unknown.

Preserve inequality qualifiers and unknown measurement types. Do not convert IC50 into Ki or Kd, average measurements, or rank unlike assay types together. Leave units null as specified above. Retain the raw response for the lifetime of the result object so parsing and interpretation remain auditable.

Save a full readable Markdown report with accession, mode, source, UTC retrieval time, unique compound count, measurement count, and all measurement rows. Print a compact summary to stdout with accession, source, status, unique compound and measurement counts, and output paths. Label the unique count as “Unique ligands listed in BindingDB for ACCESSION_ID: N”. Count distinct BindingDB monomer IDs across all returned rows, including when the same ID occurs under different target descriptions or assay types. This describes API-associated ligands, not independently verified binding partners. Include the five most frequent reported target descriptions and five most frequent measurement types, with measurement counts. Target counts count rows, not unique ligands; a ligand may occur under multiple target descriptions. Retain all target descriptions in saved reports. Sort labels by decreasing count and then lexicographically; state how many additional labels are omitted. Replace nonprintable characters in terminal labels with spaces and truncate labels to 100 characters. Do not print individual measurements to the terminal. List all measurement rows in stable order by compound ID, affinity type, affinity text, PMID, and DOI. Count unique ligands by BindingDB monomer ID without merging distinct measurements. Missing compound names must not trigger invented names or mandatory calls to another database.

Use the same stable ordering in the JSON record list and TSV. Compare all sort keys lexicographically as strings, treating null as an empty string; ties retain API order. Render Markdown with a heading `# LigandScout: ACCESSION_ID`, source, status, retrieval time, both counts, and columns `Compound`, `Reported target`, `Measurement`, `Value`, `Unit`, `References`. Link compounds to `source_url`; generate publication links only for a digits-only PMID and a DOI beginning with `10.`. Encode identifier path characters before constructing publication URLs; otherwise display escaped identifier text. Never use supplied identifiers as arbitrary URLs. Escape backslashes, pipes, HTML characters, and Markdown delimiters in cell text; replace CR/LF with spaces. Null references display `—`.

On source failure, render status `failed`, its concise error message, zero current counts, and the empty table. Never print the no-match sentence for a failed request. Persist `report.md` and print the compact summary to stdout only after all output files have been saved successfully. The API error is also logged to stderr and the command exits 1.

For no matches, print “No BindingDB records found for ACCESSION_ID.” Clarify that this does not establish the absence of binding partners. The broader target summary described in README remains future work.

### Evidence interpretation

Both the terminal summary and a nonempty saved report must state that results are **API-associated records, not confirmed binding partners**, and that target accession and assay evidence have not been independently verified. Reported target labels are source descriptions, not validation results. Preserve all rows without assigning a verified status or filtering on target names.

A BindingDB dataset DOI can contain measurements for several targets. Review the exact compound–target measurement and its accession/sequence mapping, organism, construct or processed protein, assay description, and primary reference. An IC50 describes an assay effect and does not alone establish physical binding; a Kd, a suitable binding assay, or structural evidence can support direct binding, subject to assay quality and target identity. Dataset DOI presence or target-name similarity alone cannot confirm the interaction.

Known review case: [entry 912](https://www.bindingdb.org/entry/912) contains cathepsin L rows and P0DTC1 rows. The [BDBM298298 measurement](https://www.bindingdb.org/rwd/jsp/dbsearch/Summary_ki.jsp?entryid=912&ki_result_id=352872&reactant_set_id=352874) lists P0DTC1, a 3CL-protease inhibitory assay, and IC50 2400 nM, while the entry-listing row links to an inconsistent JAK assay description. This is a data-quality flag requiring primary-source review, not grounds to infer that every row is valid or invalid. The current API response lacks row-level target accessions and stable measurement IDs, so automated evidence verification remains unimplemented. Do not scrape web pages automatically during searches or introduce a bulk-download fallback.

### Saved results

Save results under `results/ACCESSION_ID/`, relative to the current working directory, using the normalized accession. Do not create timestamp or version subfolders. Repeating a search for the same accession overwrites its generated results without prompting.

```text
results/P00699/
├── sources/
│   ├── bindingdb.json
│   └── bindingdb.raw.json
├── index.html
├── bindingdb-ligands.html
├── ligands.tsv
└── report.md
```

- `sources/bindingdb.json`: normalized source result with `schema_version`, accession, source, retrieval time, request URL/parameters, source version, status, records, and error details. Status is `success`, `no_matches`, or `failed`. Missing source versions and optional fields are null.
- `sources/bindingdb.raw.json`: original JSON response, preserved without normalization. If the response is empty or is not JSON, preserve its original bytes in `bindingdb.raw.txt` instead. Remove the obsolete raw-response variant on rerun; if no response was received, retain neither variant.
- `ligands.tsv`: compiled measurement table using the normalized source records and preserving source attribution. Initially it contains only BindingDB records; additional sources will contribute through the same normalized result structure.
- `report.md`: the full readable report; stdout contains only a compact summary.

Always make a fresh API request on rerun; saved results are output artifacts, not an implicit cache. Overwrite the source result and regenerate the compiled table and report from the current invocation. For no matches, write an empty records list and a header-only TSV. For a service or parsing failure, write status `failed`, an empty records list, and error details; regenerate the report and header-only TSV so previous successful records are not presented as current. Output-writing failures exit with code 1 and must not be reported as a completed search.

Create directories with `parents=True, exist_ok=True`. Render all artifacts in memory, stage temporary files in their destination directories, then publish with `os.replace`: current raw body first (if present), remove the obsolete raw variant, replace `ligands.tsv`, replace `report.md`, and replace `sources/bindingdb.json` last. Clean up remaining temporary files in `finally`. With no response body, remove both old raw variants. Replace only these generated filenames; preserve unrelated files and other sources' outputs. Do not load other sources' previous outputs into the compiled table in this increment.

Atomicity is per file, not across the directory. On an output failure, log that the folder may contain incomplete or mixed results, exit 1, and do not print a completed report. Do not claim rollback or use a previous report as a successful fallback. Concurrent searches for the same accession folder are unsupported in this increment. On interruption, clean up temporary files and exit 130; if publication had started, the same partial-output limitation applies. Retrieval timestamps remain inside result files even though the directory contains no run history. Help and invalid input create no result files.

## Network handling and exit codes

Use an injectable synchronous HTTPX client configured with `httpx.Timeout(30.0)` for connect, read, write, and pool operations. These are operation/inactivity timeouts, not a 30-second wall-clock deadline for the entire request; see [HTTPX timeout semantics](https://www.python-httpx.org/advanced/timeouts/). Close owned clients through a context manager; do not close a caller-injected client.

Retry `httpx.TimeoutException`, `httpx.NetworkError`, HTTP 429, and HTTP 500–599 up to three total attempts. Each attempt starts from the original endpoint, not the last redirect destination. Use one retry loop in the adapter; disable transport-level retries to avoid multiplying attempts. Without a valid `Retry-After`, wait 1 second before attempt 2 and 2 seconds before attempt 3, without jitter. Do not sleep after the final attempt. Ambiguous empty HTTP 200 bodies and decoded JSON empty strings also retry within the same attempt budget, honoring Retry-After. Other HTTP statuses, other schema errors, protocol errors, and redirect-policy failures are not retried. Non-retryable HTTPX request exceptions are recorded as transport failures; explicit redirect-policy failures use the HTTP classification above.

For retryable HTTP statuses, parse `Retry-After` as nonnegative integer seconds or an HTTP date relative to the injected UTC clock. A past date means zero delay; malformed values fall back to exponential backoff. Wait at least the exponential delay and at least the valid requested delay. If the resulting delay exceeds 60 seconds, stop immediately with an HTTP failure explaining that the service should be retried later. Inject the sleeper so tests never actually wait. Never log credentials or entire payloads at `INFO`.

| Exit code | Meaning |
| --- | --- |
| `0` | Search completed, including a successful search with no matches; help displayed successfully |
| `1` | Network, HTTP, response-schema, or output-writing failure |
| `2` | Invalid CLI input |
| `130` | User interrupted the operation |

## Modules

| Module | Responsibility |
| --- | --- |
| `cli.py`, `__main__.py` | Typer interface, logging setup, entry points, and exit handling |
| `accessions.py` | Accession normalization and syntax validation |
| `bindingdb.py` | HTTP requests, retry policy, response validation, and field mapping |
| `models.py` | Measurement, provenance, and search-result models |
| `search.py` | Search orchestration |
| `report.py` | Stable human-readable output |
| `storage.py` | Per-source JSON/raw responses, compiled TSV, Markdown report, and overwrite handling |

Keep HTTP access and rendering separate so future API adapters can reuse the command workflow. Do not add classifier callbacks, additional source integrations, bulk-input options, or local database search in this increment.

Use standard-library dataclasses for `Measurement`, `SourceError`, `SourceResult`, and `SearchOutcome`; no model-validation dependency is needed. `SearchOutcome` holds the source result and optional complete response bytes. Keep these interfaces independently testable:

```text
normalize_accession(value: str) -> str                 # ValueError on invalid input
fetch_bindingdb(accession, *, client, sleep, now) -> SearchOutcome
run_search(accession, *, client, sleep, now) -> SearchOutcome
render_summary(result: SourceResult, *, folder: Path) -> str
render_report(result: SourceResult) -> str             # Markdown ending in newline
render_tsv(result: SourceResult) -> str                # TSV ending in newline
save_results(outcome: SearchOutcome, *, root: Path) -> Path
```

`root` is `Path.cwd() / "results"` in the CLI and a temporary root in tests. `save_results` renders and publishes all artifacts, returning the accession directory; it raises `OSError` for output failures. `fetch_bindingdb` returns source failures as data and lets `KeyboardInterrupt` propagate. `run_search` orchestrates BindingDB only and does not write files. The CLI owns client creation, calls `run_search`, saves the outcome, prints the compact summary, and maps the result status or exception to an exit code. Do not catch unexpected programming errors as no matches.

## Required pytest coverage and acceptance

Use mocked HTTP responses and checked-in, small response fixtures for normal tests. Do not depend on live services or fixed live result counts.

- CLI: root positional syntax; equivalent short/long help; missing or invalid accession; normalization; unknown flags are rejected without side effects; log-level validation; stdout/stderr separation; both installed and module entry points.
- API: exact URL and parameters; no cutoff; positive responses and multiple measurements per compound; optional fields; unknown additional fields; original affinity qualifiers; both documented and observed empty responses.
- Failures: malformed JSON, HTML, missing/wrong envelope, invalid records, HTTP errors, timeout/connection retries, exhausted retries, and `Retry-After` handling. Inject the sleeper so tests do not wait.
- Report: bounded terminal summaries with target/type counts, full saved reports without row truncation, explicit unverified-evidence language, stable ordering, distinct compound and measurement counts, references, missing units, and a clear no-match message.
- Shared schema: explicit null fields, source status/error metadata, publication identifiers, duplicate measurement preservation, and a TSV round trip with embedded tabs/newlines. Confirm that report rendering escapes Markdown table delimiters.
- Workflow: a mocked accession-to-report success, no-match result, failure exit codes, and interrupt handling.
- Storage: use pytest temporary directories to verify the accession folder layout, normalized JSON metadata/statuses, raw-response variants, TSV/Markdown output, and rerun overwrites. Cover success followed by no matches or failure, obsolete raw-file removal, output-writing errors, and preservation of unrelated files. Verify that help and invalid input create no files.

Use `httpx.MockTransport` for response and exception scenarios, a fixed UTC clock, and a list-collecting sleeper. Fixtures must cover the positive envelope above, an empty list, an empty body, a JSON empty string, invalid required/optional field types, and HTML with HTTP 200. Assert the complete request parameters, saved schema keys/types, and deterministic report/TSV rows. Capture the warning and nonzero exit for a simulated failure during publication; do not test only the happy path.

Additional boundary checks: six- and ten-character accessions; whitespace and case normalization; no-match success after an earlier successful run with records; HTTP 404 with an empty body remains a failure; a last-attempt transport failure removes any old raw body; a supported HTTPS redirect preserves source query provenance; insecure/malformed redirects fail; a redirect loop stops after five hops; the retry counter is not reset by redirects; unknown measurement types survive; and a failed source never emits a no-match message.

Implementation sequence:

1. Add the HTTPX dependency, models, accession validation, and parser tests.
2. Add the API client and deterministic request/retry tests.
3. Add report/TSV rendering and overwrite/storage tests.
4. Connect the CLI to search orchestration and add end-to-end mocked CLI tests. Retain no-argument help behavior.
5. Update README with the command, result layout, overwrite behavior, and current source coverage. Run `python -m pytest -q` after an editable install with test dependencies, plus `git diff --check`.

Run the relevant pytest suite when implementing or changing code and report its results. Optional live smoke checks may use P00533 as a previously observed positive case and P00699 as a previously observed no-match case; database contents may change. P00176 returned an empty response during implementation verification, so it is not a reliable positive control. Empty unstructured responses, including those for P00699, are indeterminate rather than confirmed no-match cases. Completion requires working packaging/entry points, passing tests, and README command examples consistent with this specification.

### Local HTML ligand reports and database summary

Each search also writes `bindingdb-ligands.html` beside `report.md`, staged and
published using the same file lifecycle as the other outputs. The Markdown report
contains a `Database | Ligands Count | Ligands Report` table with a relative link
to this HTML file. There is currently one row, BindingDB; other database adapters
are not yet implemented. Counts are distinct source compound IDs, not measurement
counts. Failed searches display `Unavailable` in the table, rather than suggesting
a successful zero-match search.

The standalone HTML requires no network assets or JavaScript. Group records by
compound ID, sort IDs lexicographically, and render one row per ligand with supplied
names (or `Not supplied`), a visible SMILES column with each distinct nonblank
supplied string in a separate selectable, wrapping monospace block (or `Not supplied`),
distinct reported target descriptions, a compound link,
and expandable measurement details containing every assay value, target, and
available PMID/DOI. Escape all remote text and build BindingDB compound URLs from
encoded IDs. Include accession, source, retrieval time, status, unique and measurement
counts, evidence limitations, and links to the Markdown and TSV. Empty and failed
searches replace previous HTML with the current status and no old ligand rows.

Target identity remains unverified: the current `getLigandsByUniprots` response
supplies descriptive `query` labels but no per-measurement target accession or
sequence. Do not discard records by comparing descriptions to protein names.
Reliable future filtering requires assay-linked target accessions or sequences
resolved through online APIs, explicit treatment of fragments and mutants, and
separate matched, mismatched, and unverified records with exclusion provenance.
