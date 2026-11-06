# Completed query reuse and `--force`

The CLI accepts `--force`, a boolean flag defaulting to false. Inputs are the normalized Swiss-Prot accession and `results/` under the current working directory. Queries continue to use online APIs; this feature reuses saved query outputs rather than a downloaded database. Library callers of `run_search` may supply `root` to enable reuse; omitting it retains fresh-query behavior.

For BindingDB, UniProt and RCSB PDB independently, before querying:

1. With `--force`, bypass all saved-result checks and execute the usual API workflow.
2. Otherwise, read `results/ACCESSION/sources/SOURCE.json` with the existing versioned result reader. Require the requested source/accession, a supported schema and selection policy, status `success` or `no_matches`, no error, and complete counts. Failed or incomplete searches are never reused.
3. Require the corresponding nonempty `SOURCE-ligands.html`. For BindingDB, require exactly one saved raw response variant and verify that parsing it reproduces the normalized records. For UniProt and RCSB, require the matching accession request manifest, numbered exchanges, valid request IDs, and all response bodies with matching SHA-256 hashes. Manifest raw paths must follow the owned numbered-file layout and must not be symlinks. Any recorded request references in normalized results must exist in the manifest.
4. Missing, malformed, incompatible or inconsistent artifacts cause that source alone to be queried again. Log the reason at INFO. A valid saved source is skipped, with an INFO message containing its original retrieval date. No automatic expiry applies; use `--force` to refresh older results or successful results with optional enrichment warnings.

Run sources in the existing order. Reused UniProt context remains available to a fresh RCSB query. Preserve reused source JSON, raw files and manifests verbatim, including their timestamps; regenerate source HTML with the current renderer while retaining the original retrieval date; no provenance is fabricated for skipped requests. Mark reused outcomes in the terminal summary with the saved retrieval date. This runtime marker is not part of the saved schema.

Always rebuild all source HTML pages, `index.html`, `report.md` and `ligands.tsv` from all current outcomes, including reused sources. Missing combined reports alone do not trigger API queries. Refreshed sources retain existing publication, raw cleanup, error and exit-code behavior. `--force` can replace a previously successful snapshot with a failed or empty search, as ordinary fresh searches do. There is no directory-wide transaction or support for concurrent writers to one accession folder.

Regression coverage includes default reuse without requests, successful zero-match reuse, forced refresh, independent source invalidation, malformed data, accession mismatch, manifest corruption, missing raw/source reports, and preservation of source snapshots while rebuilding combined reports.
