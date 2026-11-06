# GO annotations in the index

## Inputs and scope

For each reviewed UniProt entry resolved from the requested accession, read all `GO` entries in `uniProtKBCrossReferences`. Use each `GoTerm` property for its name and aspect (`F` → MF, `P` → BP, `C` → CC). Deduplicate by GO ID. This section describes protein function and location; it does not change ligand selection, counts, or structural evidence.

UniProt's compact GO cross-references omit gene-product-to-term relations. Retrieve explicit relations through the [QuickGO annotation API](https://www.ebi.ac.uk/QuickGO/api/index.html), using the resolved primary accession as `geneProductId=UniProtKB:ACCESSION`, `limit=200`, and successive `page` values. This is GOA enrichment for UniProt annotations, not another ligand database integration. Do not expand the displayed term set to QuickGO-only annotations.

## Algorithm and failures

1. Extract and validate UniProt GO identifiers, names and aspects. With no terms, save `no_matches` without a QuickGO request.
2. Query QuickGO with the existing audited JSON HTTP helper, retaining every response, request parameters, retrieval timestamps and available headers in the UniProt request manifest/raw files. Existing request spacing, retries, timeouts and redirect restrictions apply.
3. Validate pagination metadata, stable hit counts, exact gene-product identity and agreement of GO aspects for matching terms. Retrieve all pages without a result cap; log annotations retrieved against the reported total after each page. Require the final row count to equal `numberOfHits`.
4. For each UniProt GO term, collect distinct explicit `qualifier` strings from matching QuickGO annotations. Preserve all relations, including `NOT` modifiers; sort relations for deterministic display. Do not assume `enables`, and do not replace other relations with it.
5. A failed or malformed optional GO lookup does not fail the ligand query. Retain the UniProt terms, mark relation lookup unavailable, and discard partially retrieved relation assignments. Distinguish this from a successful empty UniProt GO list. No evidence-level or taxon filtering beyond exact gene-product identity is applied; multiple annotations for a term are summarized as a set of relations, not a new biological assertion.

## Persistence and presentation

Add optional `source_context.go_annotations` to the UniProt saved result, with `status` (`success`, `no_matches`, `unavailable`), `terms`, `request_ids`, and nullable `error`. Each term has `id`, `name`, `aspect` and `relations`. Validate this additive context on read; retain compatibility with older results. Existing entry provenance establishes the source of the UniProt GO list; `request_ids` refer to QuickGO calls in the same UniProt manifest. QuickGO release information is retained only when supplied in captured response headers; no release is invented.

Append a GO annotations section to `index.html`, grouped into Molecular Function (MF), Biological Process (BP), and Cellular Component (CC), always in that order. Render one row per GO term in a two-column table: Description (GO identifier link and name) and Relation (explicit labels). Use the same 72%/28% column widths and fixed table layout in all three aspect tables to align the columns. Escape all text. Display `not supplied` where no explicit relation is available and explain lookup failures. Show the original source retrieval date and attribution to UniProt and QuickGO/GOA. An unavailable UniProt source must show unavailable rather than an empty annotation claim.

Completed UniProt snapshots that predate this field require a UniProt refresh on the next ordinary invocation; other completed sources remain reusable. Subsequent normal runs reuse saved GO data and rebuild the HTML offline. A successful UniProt result with unavailable optional GO enrichment is reusable; `--force` retries enrichment.

Tests cover MF/BP/CC grouping, one line per term, multiple and negated relations, absent qualifiers, pagination, exact identity/aspect validation, incomplete/service failures, empty lists, HTML escaping, saved schema validation, request manifests, legacy cache invalidation and reuse. A live smoke check for P00519 successfully retrieved 110 UniProt GO terms (27 MF, 66 BP, 17 CC) and two QuickGO pages; these counts are observations, not fixed expectations.

### Entry review status in the heading

Derive the GO heading from the saved UniProt `source_context.entry.entry_type`, originally obtained directly from the accession API's `entryType`. Show `GO annotations (UniProt | SwissProt)` for `UniProtKB reviewed (Swiss-Prot)` and `GO annotations (UniProt | <s>SwissProt</s>)` for an explicitly `UniProtKB unreviewed (TrEMBL)` entry. Use semantic HTML strike-through with an unreviewed-entry tooltip. Absent/unrecognized metadata or a missing source shows `SwissProt status unavailable`; do not equate unavailable with unreviewed. Clarify that this describes entry review status, not the evidence level of every GO annotation. The current adapter accepts reviewed entries only; the unreviewed rendering branch does not expand query eligibility. Direct UniProt metadata is available for all accepted entries, so no `profunc browse accession` subprocess is needed. Rendering remains offline, including when saved results are reused.
