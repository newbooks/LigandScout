# UniProt ligand annotations and binding sites

Status: implemented. Add UniProt alongside BindingDB in the accession search workflow. This retrieves annotations without calling a classifier or prediction service. Preserve the existing command, BindingDB selection, and API-only architecture.

Current shared serialization is schema v5, as defined in the [RCSB integration spec](rcsb-integration.md). The UniProt adapter and its annotation validator retain the v4 contract below internally; shared upgrades add null structural counters without changing annotations. Readers accept UniProt v4/v5.

## Scope and evidence

Query only the supplied accession's reviewed UniProt entry. Collect explicitly identified binding-site ligands and cofactors; retain generic binding sites, catalytic reactions, function, and subunit information as context. Do not assign every reaction participant as a binding partner. Do not infer ligands from names, sequence similarity, GO terms, keywords, domains, or catalytic residues.

Automatically attempt structure enrichment for each eligible identifier; “optional” below means that failure to enrich does not fail the annotation search, not that a new CLI flag is required. Use only the selected resources: UniProt for annotations and RCSB PDB's Chemical Component Dictionary (CCD) for optional SMILES enrichment. The CCD call supplies chemical identity, not evidence that the molecule binds this protein. No ChEBI, Rhea, PubChem, ChEMBL, or other external API calls are introduced; their identifiers/links may be preserved as UniProt provenance. Experimental PDB ligand contacts are retrieved by the separate implemented [RCSB adapter](rcsb-integration.md).

Present the result as “UniProt-annotated ligands and cofactors,” with annotation evidence attached. Reviewed status does not imply that every annotation was measured experimentally. Preserve ECO codes and their source identifiers; distinguish experimental, similarity-based, and other annotations without assigning invented confidence scores. No measured affinity values are fabricated.

## Requests and entry identity

```text
GET https://rest.uniprot.org/uniprotkb/ACCESSION.json
Accept: application/json
User-Agent: ligandscout/VERSION
```

Retrieve the complete entry, including sequence, `features`, `comments`, references, and entry audit metadata. Retain `primaryAccession`, `secondaryAccessions`, `entryType`, `uniProtkbId`, protein description, organism/taxonomy, sequence length and sequence, entry/sequence versions and last update dates. Record release headers when available; absence is not a failure.

- Primary accession equal to the query: accept.
- Different primary accession with query explicitly in `secondaryAccessions`: accept as a resolved secondary accession, preserve both IDs, and display the resolution. Keep the original normalized query as the output folder name and request provenance; use the primary accession as the annotated target ID.
- Otherwise mismatched identity: source failure. Do not silently substitute a related protein.
- Require reviewed `entryType` (`UniProtKB reviewed (Swiss-Prot)`). An unreviewed result produces an unsupported-entry source failure, not zero ligands.
- HTTP 404: source failure explaining that the accession was not found. Empty/whitespace HTTP 200 bodies and decoded JSON empty strings are retried, then fail as indeterminate; malformed nonempty JSON fails without retry.
- Missing `features` or `comments` means an empty collection; present collections must be lists with valid objects. Missing required entry identity or sequence data is a schema failure.

Use the shared synchronous HTTPX client with TLS verification and 30-second operation timeouts. Reuse BindingDB's three-attempt retry/backoff/Retry-After and HTTPS redirect rules, implemented in `api_http.py` for the new requests. Reuse existing retry-delay logic where practical, but leave the BindingDB request/parser implementation unchanged in this increment; its tests remain the compatibility contract. Retry transient transport failures, HTTP 429/5xx, and ambiguous empty bodies/decoded JSON empty strings; do not retry 404 or malformed nonempty responses. Follow only explicit HTTPS redirects, up to five per attempt, with retries restarting the original URL.

UniProt entry retrieval is one object, so it has no pagination. BindingDB, UniProt and RCSB run independently: a failure of one source must not prevent output from the others. The CLI exits 1 if any core source fails, otherwise 0. Optional chemical-enrichment failures are warnings and do not invalidate the core UniProt result.

## Annotation extraction

| Input | Action |
| --- | --- |
| `features[].type == "Binding site"` with `ligand` | Extract ligand identity, location, ligand part, note, feature ID, sequence identifier, and evidence |
| `comments[].commentType == "COFACTOR"` | Extract each named cofactor and its ChEBI reference, evidence, and comment notes; no site is invented |
| Binding site with generic `substrate`, missing ligand identity, or unspecified ligand | Save in a separate generic-site context section; do not create an identified molecule |
| `CATALYTIC ACTIVITY` comments | Preserve full reaction text, EC number, cross-references, and evidence as reaction context; do not turn all ChEBI cross-references into ligand rows |
| `FUNCTION`, `SUBUNIT`, `ACTIVITY REGULATION`, `BIOPHYSICOCHEMICAL PROPERTIES` | Preserve source text and evidence in optional context sections, with no free-text molecule extraction |
| `Chain` and `Peptide` features | Retain processed-product names and intervals to interpret sites on polyproteins |
| `Active site`, disulfides, domains, other features | Do not interpret as ligand-binding sites |

Use case-sensitive feature/comment types from the JSON contract. Ignore unknown feature/comment types and additional fields for forward compatibility. An absent evidence list means no evidence supplied, not experimental support. An absent optional ligand-part, feature-ID, label, or note field remains null. Malformed recognized ligand/site/cofactor annotations fail the core source rather than silently omitting evidence. Malformed optional context produces a warning with its JSON path and retains the raw entry.

### Identity and grouping

Normalize `ChEBI:CHEBI:30616` and `CHEBI:30616` to `CHEBI:30616`. Validate positive numeric ChEBI IDs; retain the original ID and JSON path. Use the typed `ligand.id` as the whole-ligand identifier. Feature cross-references may also describe a ligand part; retain them under their original annotation and apply the explicit rules below rather than requiring every ChEBI reference to equal the whole ligand.

Group the same ChEBI ID across features and cofactor comments into one record containing all roles, sites, and evidence. Do not merge different ChEBI IDs, protonation states, stereoisomers, or ligand parts. Preserve separate binding events and `ligand.label` values; grouping does not imply one physical binding pocket.

For explicit names with no chemical ID, assign a reproducible local ID `UNIPROT:PRIMARY:NAME:SHA256`, using the full SHA-256 hex digest of the casefolded, whitespace-collapsed name. Merge only identical normalized names in this ID-less group. Do not merge ID-less names with identified compounds by spelling alone. Names `substrate`, `ligand`, `unknown`, or `unspecified` (case-insensitive) are generic context, not molecule IDs. Names associated with a valid ChEBI ID remain identified even when broad chemical classes; label such records as annotations and do not presume a defined structure.

### Binding-site coordinates

- Preserve full `location.start`/`end` objects, including `value` and `modifier`, and `location.sequence` when supplied. Coordinates are **1-based, inclusive UniProt sequence coordinates**, never PDB residue numbers.
- An exact single residue can show `271 (K)`; an exact interval shows `248–256` with its sequence segment in the details. Derive residues only when the site refers to the retrieved canonical sequence and positions are in bounds. Verify against that sequence.
- Unknown, uncertain, before/after, or other non-exact endpoints retain their modifiers and display as uncertain; do not turn them into exact residues. Invalid exact intervals (start greater than end, zero, or beyond sequence length) are schema failures.
- A site on an explicitly identified alternate sequence is retained with that identifier; do not map it to canonical residues or retrieve isoforms automatically. Show “alternate sequence; residue letters unavailable.”
- For canonical exact sites, list every containing processed `Chain`/`Peptide` feature whose exact coordinates contain the whole interval. Keep parent coordinates; do not select a single product when none matches or overlapping products make the assignment ambiguous.
- Preserve `ligandPart` independently (e.g. the bound portion of ATP) and feature evidence. A phosphate-binding part does not identify phosphate as a separate full ligand.
- Cofactor-only records display “Binding site not annotated.” Generic sites appear separately from ligand-linked sites. Do not connect sites to cofactors solely because names seem similar.

## SMILES enrichment using RCSB CCD

UniProt may supply ligand identifiers without SMILES. For each unique normalized ChEBI ID, search CCD by the exact external accession; ID-less or generic annotations skip enrichment.

Verified request example (November 5, 2026):

```text
GET https://search.rcsb.org/rcsbsearch/v2/query?json=URL_ENCODED_JSON
```

```json
{
  "query": {
    "type": "terminal",
    "service": "text_chem",
    "parameters": {
      "attribute": "rcsb_chem_comp_related.resource_accession_code",
      "operator": "exact_match",
      "value": "CHEBI:15422"
    }
  },
  "return_type": "mol_definition",
  "request_options": {"paginate": {"start": 0, "rows": 100}}
}
```

Paginate until the validated `total_count` is exhausted; inconsistent totals, repeated IDs, or incomplete pages fail that ligand's enrichment with a warning. HTTP 204 is a documented search no-hit outcome; HTTP 200 must have a valid structured search response. For every candidate ID, request `GET https://data.rcsb.org/rest/v1/core/chemcomp/CCD_ID` and verify its identity and `rcsb_chem_comp_related` contains the pair `resource_name == "ChEBI"` and `resource_accession_code == requested ChEBI ID` in the **same object**. A match from a different resource is not acceptable.

Read `rcsb_chem_comp_descriptor.SMILES_stereo` preferentially, falling back to `SMILES` with an explicit non-stereospecific label. Read `InChIKey` if supplied. Preserve strings verbatim and record CCD ID, descriptor field, retrieval time, request ID, and the exact ChEBI cross-reference mapping. This is cross-reference-based chemical enrichment, not an experimentally validated ligand assignment.

- One distinct eligible descriptor: populate `smiles` and its provenance only after all search pages and candidate records have been retrieved and validated.
- Several CCD candidates with identical selected descriptor text, descriptor field, and InChIKey (including all-null InChIKeys): use that descriptor and preserve all candidate IDs.
- Different selected descriptor text, descriptor fields, or InChIKeys: leave the main SMILES null and display candidate structures as ambiguous alternatives with links. Do not select alphabetically.
- No mapping, absent structure, ID-less annotation, or unavailable service: leave SMILES null with a reason (`no_exact_mapping`, `no_descriptor`, `ambiguous`, `unavailable`, or `not_applicable`). A service failure must not masquerade as no mapping.
- Missing sites do not prevent structure enrichment, and missing structures do not prevent site reporting.

Do not search by name, replace a ligand with its parent/neutral form, build SMILES from a formula, predict stereochemistry, or copy a BindingDB structure into a UniProt row. Do not follow reaction cross-references to enumerate additional compounds.

Important verified example: P00519 annotates ATP as `CHEBI:30616`, whereas RCSB component `ATP` cross-references `CHEBI:15422`. These IDs are different; the ATP name or CCD code must not be used to bridge them. An exact search may find another valid component or yield no structure, which must be reported honestly.

If any candidate/page request fails or returns an invalid response, mark the whole chemical lookup `unavailable`, retain any alternatives already observed, and do not select a partial result. A candidate fetched successfully but without the exact paired ChEBI mapping is discarded, with its reason retained. All candidates discarded means `no_exact_mapping`. An exact-mapped candidate without a descriptor prevents selecting another candidate as uniquely resolved: use `no_descriptor`, retaining any alternatives.

Cache CCD enrichment by ChEBI ID only within the current run (including no-hit and failure outcomes). Use sequential requests at no more than two starts per second across all UniProt/CCD attempts and redirect hops, an injectable monotonic clock and sleeper, and the core retry rules. Do not make a top-level RCSB report row for these chemical lookups: the interaction source remains UniProt, and RCSB is identified as the structure provider.

## Shared models, counts, and persistence

Introduce saved schema **version 4**, while reading BindingDB v1/v2/v3 and retaining its record ordering, selection, and existing TSV column prefix. Add these optional fields to `Measurement` with empty/default values for legacy BindingDB records:

- `ligand_roles: list[str]`: `binding_site` and/or `cofactor`.
- `binding_sites: list[dict]`: location, sequence identifier, residue/segment if derivable, ligand part/label, description, feature ID, JSON path, processed-product context, and evidence.
- `annotation_evidence: list[dict]`: annotation origin/path, ECO code, source, and source ID; preserve association with the site or cofactor comment instead of flattening evidence into a confidence label.
- `structure_provenance: dict`: enrichment status/reason, requested chemical ID, chosen descriptor field, candidate CCD records and request references. Empty for old records.

Add `source_context: dict` to `SourceResult` (default empty): entry identity/version, sequence context, generic sites, reaction descriptions, function/subunit/regulation/biophysical comments, and extraction counters. BindingDB gets the empty default and upgrades to v4 without invented UniProt metadata. Add source `UniProt` and policy `uniprot-annotated-ligands-v1` to validated source policies. Extend the saved-result validator for the new fields and source-specific invariants; do not route v1 data through a UniProt upgrader.

One UniProt record represents one normalized ligand group with all its annotations, **not an affinity measurement**. Populate `compound_id`/`compound_name` from ligand identity and `target_id`/`target_accessions` from the resolved protein entry, retain organism/taxonomy, `source_url` as the UniProt entry link, `evidence_type="uniprot_annotation"`, and nullable `smiles`/`inchikey`. Leave affinity values, units, assay IDs, target confidence, and mutation fields null. Set `target_match_status="entry_annotation"`; do not claim direct assay assignment. Individual PubMed IDs stay in annotation evidence; never pick an arbitrary first paper as the ligand's sole reference.

UniProt counts retain the shared v4 count fields, with activity-specific fields null. Add source-specific fields `record_kind="ligand_annotation"`, `annotation_records`, `unique_ligands`, `binding_site_annotations`, `generic_binding_sites`, `cofactor_annotations`, and structure-availability counts. Counts are nonnegative on a complete source; complete-result counts derived from records must agree. Failed-source complete counts are null. Change the reader's activity-count-equals-record-count invariant to apply to BindingDB only. The unique-ligand count includes explicit name-only groups and separately identifies how many have chemical IDs.

- Identified records present: `success`.
- Valid entry with no identified ligand groups: `no_matches`, meaning “No specifically identified ligands/cofactors annotated in UniProt.” Still show any generic sites and functional/reaction context; this is not proof of no ligands.
- Core retrieval/extraction failure: `failed`; do not expose partial ligand records as complete. Preserve raw responses and error details.
- CCD failures: core status remains complete with warnings, unknown structures, and explicit reasons.

Save alongside BindingDB under `results/QUERY/`:

```text
sources/uniprot.json
sources/uniprot-requests.json
sources/uniprot-raw/000001.json
uniprot-ligands.html
index.html
report.md
ligands.tsv
```

Store each complete UniProt/CCD HTTP response separately (`.txt` for non-JSON). The request manifest records monotonically increasing exchange ID, role, query chemical IDs, URL/parameters, attempt/redirect, UTC time, HTTP status or transport error, selected response headers, relative raw path, and SHA-256. Record empty final responses and enrichment failures. Retain core entry JSON, not only extracted fields. Never fabricate a database release from the retrieval date.

Stage artifacts before publishing through the current writer; overwrite current generated source files, source HTML, index, and combined reports on rerun. Remove only stale numbered raw files owned by this manifest after successful publication, preserving unrelated files and historical removed-source artifacts. No timestamp subdirectory. Completed source snapshots are reused by default under the shared [query reuse policy](query-reuse.md); `--force` bypasses reuse. Only refreshed source artifacts are overwritten. Keep per-file atomicity and existing output failure/interruption behavior.

## Presentation beside BindingDB

`run_search` returns `[BindingDB, UniProt]` in that order even if one fails. Index and Markdown summaries contain one row per queried source, separate ligand counts, status and retrieval times, and links to each source's local HTML and JSON. Never sum cross-source ligand counts or link to an unimplemented PDB interaction search.

The summary/index “Measurements” heading becomes “Records,” with source-specific explanatory labels: BindingDB measurement records versus UniProt ligand-annotation records. BindingDB terminal details may keep “Measurements”; UniProt uses “Ligand annotation records.” BindingDB detail tables keep their measurement semantics. The compact terminal output lists each source and remains bounded; it links to `index.html`, not individual HTML filenames.

The UniProt HTML adapts to available data:

1. Target identity, reviewed status, resolved accession, organism, sequence length, entry/sequence versions.
2. An annotated-ligand table when records exist: name, chemical ID/link, role(s), selectable SMILES or explicit missing reason, structure provider, and binding-site summary.
3. Binding sites displayed as bullet lists with coordinates/residues, sequence/processed-product context, and ligand parts. Detailed annotation evidence remains in the linked JSON; cofactor prose remains in context sections and CCD information in the structure cell.
4. Separate generic-site and reaction-context sections when available; function, subunit and regulation sections only when present. If only context exists, show it with the no-identified-ligands message.
5. Warnings for uncertain coordinates, alternate sequences, inferred annotations, missing structures, and source/enrichment failures. Do not label annotation evidence “confirmed binding.”

Escape all remote text in HTML/Markdown, validate/encode external identifiers before building links, and use relative links back to the index. Render locally without scripts or third-party assets. Append the new annotation/site/structure fields to TSV as compact JSON; leave measurement columns blank for UniProt rows. Do not infer shared identities or transfer sites across sources. Clearly label `smiles` as supplied by CCD and show chemical-match uncertainty.

### Human-readable HTML and source links

The UniProt HTML must not embed JSON objects, raw annotation paths, or `<pre>` data dumps. Ligand rows have four columns: Ligand, Role, SMILES / structure provider, and Binding sites. Do not include an Annotation details column. Render identified and generic binding sites as HTML bullet lists, one item per site; use “Binding site not annotated” when absent. Keep readable structure-status/CCD links in the structure cell, cofactor prose in its context section, and detailed evidence in the linked saved JSON. The complete saved JSON and request manifest remain available as linked provenance.

“Function context” comes from UniProt `comments` with `commentType="FUNCTION"`, specifically `texts[].value`. Other context kinds use their own corresponding comments. Group comments by kind in source order and generate a deterministic extract of at most three leading source sentences per section. Normalize whitespace, preserve source wording/negation and parenthetical evidence qualifiers, and avoid splitting common abbreviations. Do not use a generative service or imply that this extract covers the entire source. When more sentences exist, state that additional information is available in UniProt. Cofactor notes use the same short-extract rule. Known structured reaction fields can be rendered as labeled prose; numerical biophysical values require explicit source units.

Link the protein header, ligand table, and each context section to the resolved primary UniProt entry. Include a direct link to the UniProt JSON API record in the protein header. These live links may show a newer entry version; the local JSON remains the archived evidence for this search. Escape all source prose and notes. Regression tests must verify the four-column table and binding-site bullet lists, absence of inline JSON/pre blocks, source links, three-sentence extraction, retained qualifiers, supported structured prose, and HTML escaping.

## Exact validation and serialization contract

### Core entry and source metadata

Use the existing accession normalizer for query, primary, and secondary IDs. Require a JSON object, a valid `primaryAccession`, the reviewed entry type, a `sequence` object containing a nonempty uppercase ASCII-letter `value`, and a positive integer `length` equal to `len(value)` (reject bools as integers). Do not restrict to the 20 standard amino acids. Require `organism.scientificName` as a nonempty string and `organism.taxonId` as a positive integer. `secondaryAccessions` defaults to `[]`; a present value must be a list of accession strings. Require feature/comment items to have string `type`/`commentType` respectively, even for unknown kinds.

`proteinDescription.recommendedName.fullName.value` is the preferred display name; otherwise use `uniProtkbId`, otherwise the primary accession, with a warning. Entry-audit versions/dates, reference lists, optional organism names, and `uniProtkbId` may be absent and remain null/empty; preserve malformed optional metadata as a warning rather than inventing a value. `source_version` is the final core response's `x-uniprot-release` header or null. Preserve `entryAudit.entryVersion` and `sequenceVersion` separately from release. Use the core response completion time for `retrieved_at`; later CCD calls have their own timestamps.

The root `request` is `{"url":"https://rest.uniprot.org/uniprotkb/QUERY.json","parameters":{}}`. Core failures use existing error kinds `transport`, `http`, or `schema`; unsupported/unreviewed or mismatched entries use `schema`. `error.http_status` is the final core status if a complete response was received, otherwise null. Enrichment failures are warnings with `{code, message, role}`; the role is a chemical ID such as `CHEBI:30616`. Never embed whole response bodies in errors.

### Annotation objects

All paths use JSON Pointer syntax (`/features/0`, `/comments/2/cofactors/0`). Extract features in entry order, followed by cofactor comments/items in entry order. Retain this order inside a ligand; sort final ligand records by `compound_id` lexicographically. Preserve all source names as an encounter-ordered unique alias list in `source_values.aliases`; choose the first nonempty explicit name for display, or null if only an ID is available. Sort `ligand_roles` in the fixed order `binding_site`, `cofactor`. Name-only IDs hash UTF-8 bytes of the normalized name, with no newline.

Validate a present `ligand`/`ligandPart` as an object. `name`, `id`, `label`, `note` are optional strings or null; preserve additional fields in the raw response. Treat absent ligand objects and valid empty objects as generic annotations. Invalid ChEBI-prefixed identifiers fail the source. An explicit ID from a currently unsupported namespace is retained in `source_values` with a warning, and the explicit name becomes a name-only group; without a name the annotation stays generic. A ChEBI ID may identify a ligand even if its name is absent.

Cofactor identities come from `name` and `cofactorCrossReference` (`database`, `id`); comment notes are `comments[i].note.texts[]`. Save each original cofactor object and note with its path in `source_values.cofactor_annotations`. Do not infer stoichiometry from note text. A missing `cofactors` list means no named cofactors, with notes retained as context; a present non-list or malformed cofactor is a schema failure. Generic cofactor names stay in `source_context.generic_cofactors`, not ligand rows.

For binding-feature ChEBI cross-references, distinguish the whole-ligand ID and `ligandPart.id` explicitly. Additional different ChEBI IDs cause a warning (`unassigned_feature_cross_reference`) rather than being promoted to ligands or automatically treated as contradictory whole-ligand assertions. If `ligand.id` is absent, do not choose a whole-ligand ID from an untyped cross-reference list; use the explicit name or generic context. This avoids confusing heme with its bound iron atom.

Each evidence item is `{annotation_path, evidence_code, source, source_id, category}`. Require `evidenceCode` matching `ECO:[0-9]{7}`; optional `source` and `id` must be strings or null. Map only `ECO:0000269` to `experimental`, `ECO:0000250` to `similarity`, and `ECO:0000255` to `sequence_analysis`; all other codes are `other`, preserved verbatim. Missing evidence is an empty list labeled “Evidence not supplied,” never experimental. Do not fetch an ECO ontology or assume an unknown code's descendants. If evidence is malformed on a selected ligand/site/cofactor, fail extraction. For each feature use its `evidences`; for each cofactor concatenate its own and parent comment evidence with their original paths. Do not attach unrelated reaction evidence or all entry references to a ligand.

Each `binding_sites[]` object has exactly these keys:

```text
annotation_path: str
feature_id: str|null
location: {start: {value: int|null, modifier: str},
           end: {value: int|null, modifier: str}, sequence: str|null}
sequence_scope: "canonical"|"alternate"
residue_sequence: str|null
description: str|null
ligand_part: {name: str|null, id: str|null}|null
ligand_label: str|null
processed_products: list[{annotation_path: str, name: str|null, start: int, end: int}]
evidence: list[evidence item]
```

Keep original location/ligand-part objects under the annotation's `source_values` entry if normalization drops extra fields. Missing endpoints normalize to `{value:null, modifier:"UNKNOWN"}`. A present modifier must be a nonempty string; a positive integer value with no modifier defaults to `EXACT`; no value defaults to `UNKNOWN`. An explicit `EXACT` endpoint without a positive integer value is invalid. Unsupported modifiers are retained with a warning and treated as uncertain. Do not validate alternate-sequence coordinates against canonical length. Sequence is canonical when absent or exactly the primary accession; all other sequence identifiers are alternate. A site spanning two products belongs to neither unless one annotated product contains the entire interval. `processed_products` preserves every valid containing feature; overlapping products are not resolved arbitrarily.

For context comments, retain `{comment_type, annotation_path, data}` where `data` is the original JSON object. Render `texts[].value` and evidence for text comments; render the reaction's `name`, `ecNumber`, cross-references and evidence for catalytic activity. For biophysical properties render available source text and known numeric properties only when units are explicit; omit unsupported shapes from the human-readable HTML rather than dumping JSON or inventing units. Never traverse these fields to mine ligand names. Raw structured data remains available even if the renderer cannot summarize a new property.

### Schema-v4 shapes and invariants

All existing `Measurement` fields remain. New defaults are `ligand_roles=[]`, `binding_sites=[]`, `annotation_evidence=[]`, `structure_provenance={}`; `SourceResult.source_context={}`. Append these four measurement fields to the current TSV columns, preserving the entire existing prefix. `source_values` holds `aliases`, `binding_annotations` (original feature plus path), and `cofactor_annotations` for UniProt; it does not replace typed fields.

For UniProt set `measurement_type`, `measurement_value`, `measurement_unit`, `measurement_relation`, `measurement_upper_value`, `measurement_text`, `measurement_basis`, `pmid`, `doi`, `assay_id`, `assay_type`, `source_record_id`, `document_id`, `target_confidence`, `target_relationship`, `variant_accession`, `variant_mutation`, and `variant_isoform` to null. Set `target_description` to the selected protein display name, `target_id` to the primary accession, `target_accessions=[primary]`, and organism fields to the entry organism (taxon ID serialized as string to match the shared model). `compound_id` remains the ligand's chemical/local ID, never the protein accession. Build `source_url` from the validated primary accession only.

`source_context` has these keys for every completed UniProt result: `entry` (the normalized entry metadata above), `generic_sites` (site objects), `generic_cofactors` (original objects and paths), `context_comments` (the context objects above), and `processed_products` (validated Chain/Peptide context). On core failure use `{}` and `records=[]`. Optional malformed context objects are omitted from normalized context with a warning and preserved in raw JSON. Do not retain a partially extracted complete-looking context on failure.

Every v4 `counts` object retains the existing v3 fields plus these exact additional keys:

```text
record_kind                 # "measurement" (BindingDB) or "ligand_annotation" (UniProt)
annotation_records
identified_chemical_ids
binding_site_annotations
generic_binding_sites
cofactor_annotations
smiles_available
smiles_unavailable
```

For completed UniProt results: `complete=true`, `unique_ligands=annotation_records=len(records)`, `identified_chemical_ids` counts ChEBI-keyed records, `binding_site_annotations=sum(len(r.binding_sites))`, `generic_binding_sites=len(source_context.generic_sites)`, and `cofactor_annotations` counts extracted cofactor items including retained generic items. `smiles_available` counts records with a selected nonempty SMILES, and `smiles_unavailable=annotation_records-smiles_available`. Existing `candidate_targets`, `eligible_targets`, `candidate_activities`, `included_activities`, `excluded_activities`, `unverified_activities`, `candidate_assays`, and `eligible_assays` are all null. Failed UniProt results have `complete=false`; every numeric count, including `unique_ligands`, is null. `record_kind` remains `ligand_annotation`.

For BindingDB v4 keep all v3 values/invariants, add `record_kind="measurement"`, and initialize all seven new numeric keys to null (the existing `unique_ligands` remains its distinct compound-ID count). The reader must not run the old all-count-values-are-integers check on `record_kind`. Unknown/missing v4 count keys are validation errors; reject bool values used as counts. Keep the existing BindingDB selection audit. UniProt uses exactly `{outcome_reason, targets:[], activities:[], assays:[]}`, with reason `incomplete` on failure, `no_identified_annotations` on no matches, and null on success. Version upgrades add defaults only; they must not invent site or structure provenance.

### Structure provenance and request manifests

`structure_provenance` has `{provider, status, requested_chemical_id, smiles_field, candidate_components, request_ids}`. `provider="RCSB CCD"`; status is one of `available`, `no_exact_mapping`, `no_descriptor`, `ambiguous`, `unavailable`, `not_applicable`. `smiles_field` is null unless a single descriptor choice is established. Each candidate records CCD ID, exact matching cross-reference object, selected descriptor field/string, InChIKey, and request ID. For `available`, main `smiles` and `inchikey` must equal the selected candidate data; for all other states main `smiles`/`inchikey` are null. All ID-less records get `not_applicable` without network access. Do not infer `not_applicable` from a name when a valid chemical ID exists.

Search response `total_count` must be a nonnegative integer and each returned identifier a nonempty ASCII CCD ID (`[A-Za-z0-9]+`), URL-encoded on retrieval. Start pagination at 0 and advance by returned result count, requiring all nonfinal pages to contain 100 IDs. No duplicates across pages; terminate only when accumulated IDs equals total. A first-page 204 means no mapping; 204 after positive pages means incomplete lookup. Validate `result_type=="mol_definition"` on structured responses. Missing descriptor collections are allowed, but present malformed descriptor types fail enrichment. For a component response require `rcsb_id` and `chem_comp.id` to equal the requested CCD ID, an object `rcsb_chem_comp_descriptor` if present, and an array `rcsb_chem_comp_related` if present. Missing related mappings means no eligible exact mapping; malformed present related entries are enrichment failures. Read only the two named SMILES fields and InChIKey as strings/null; a nonempty primary descriptor wins within a component. Whitespace-only descriptor strings count as absent. Never compare different descriptors through chemical canonicalization.

`uniprot-requests.json` has exactly `{schema_version:1, accession, exchanges}`. Each exchange has `sequence` (1-based global exchange ID), `request_id` (1-based logical request ID shared by retries/hops), `role` (`entry`, `ccd_search`, `ccd_component`), `chemical_id` (null for entry), `attempt` (1-based), `redirect_hop` (0-based), `url`, `parameters`, `retrieved_at`, `http_status`, `transport_error`, `response_headers`, `raw_path`, and `sha256`. Raw paths are relative to `sources/` and use the exchange ID as `uniprot-raw/000001.json` or `.txt`. An empty received body is a real `.txt` file with a hash; no complete response means null path/hash. Selected lower-case headers are `date`, `etag`, `last-modified`, `content-type`, `retry-after`, `x-uniprot-release`, and `x-uniprot-release-date`. Never persist authorization/cookie headers. Hash the exact retained response bytes, not reserialized JSON. Add a typed HTTP-exchange dataclass and `SearchOutcome.exchanges` with default `[]`; BindingDB continues saving its final raw body and does not need a new manifest.

Stage new files, publish per-source raw bodies and HTML, TSV/Markdown, index, manifests, and normalized JSON last. Remove only obsolete `uniprot-raw/[0-9]{6,}.(json|txt)` files after successful publication; do not trust arbitrary paths from an old manifest. Keep numbered files not matching this pattern untouched. Render every artifact from current-run objects, not saved source JSON.

## Module interfaces and completion checks

```text
api_http.py:
    request_json(url, *, params, role, chemical_id, audit, client, sleep, now, monotonic)
        -> HttpResult                         # data/status or typed request failure
uniprot.py:
    parse_entry(data, *, query) -> ParsedEntry
    fetch_uniprot(query, *, client, sleep, now, monotonic) -> SearchOutcome
ccd.py:
    enrich_ligands(result, *, audit, client, sleep, now, monotonic) -> SourceResult
models.py:
    upgrade_result(result) -> SourceResult    # output v4, idempotent for v4
report.py:
    render_uniprot_html(result) -> str
search.py:
    run_search(accession, *, client, sleep, now, monotonic=time.monotonic) -> RunOutcome
```

`HttpResult` carries parsed data, complete response bytes or null, final status, selected headers, and an optional `SourceError`; a first-page CCD-search 204 is an allowed empty search result, never a UniProt-entry success. The audit collector allocates request/exchange IDs and records every attempt before errors are converted to source failures. Generic schema failures are never caught as no matches. `fetch_uniprot` owns extraction, enrichment, result timestamp/provenance, and error conversion; it does not create or close the injected client or write files. `ParsedEntry` is a dataclass with `records`, `source_context`, `counts`, `selection_audit`, and `warnings` as described above, plus `status` (`success` or `no_matches`). It contains no transport metadata. `parse_entry` raises a typed `UniProtSchemaError` for invalid core data; `fetch_uniprot` builds the complete `SourceResult` with request/retrieval metadata after parsing. This avoids constructing a partially populated saved-result object. `enrich_ligands` changes only structure fields, warnings, and structure counts, preserving annotations and core source status. Let `KeyboardInterrupt` propagate to the CLI.

Throttle each new HTTP start to `max(last_start+0.5-monotonic(),0)`, including retries and redirects. Apply retry waits before rate-limit waits; do not sum them blindly. A delay exceeding the shared 60-second Retry-After limit fails the logical request rather than blocking indefinitely. A timeout is an operation timeout, not a 30-second whole-run deadline. Core requests and enrichment share the audit/rate limiter; no new async machinery or extra runtime dependencies are required.

Source-specific report dispatch must handle `BindingDB` and `UniProt` explicitly; do not reuse BindingDB compound URLs, measurement headings, or its blanket unverified-API disclaimer for UniProt. UniProt text must say annotations may be experimental or inferred and CCD supplies chemical structure only. Both sources still caution that annotation/assay association is not by itself proof of direct binding. Include a source-specific reference renderer for annotation evidence, with only validated PubMed/UniProt IDs linked and all other evidence displayed as escaped text. Avoid rendering an unbounded list of site/evidence warnings in the terminal: show up to five warnings plus an omitted count; full lists belong in saved reports.

Minimum acceptance examples (fixtures, not fixed live-content promises):

| Case | Required result |
| --- | --- |
| P00519-like ATP at 248–256, 271, 316–322 plus Mg cofactor | Two ligand groups; three ATP sites; Mg has no invented site; ATP CCD component with a different ChEBI ID cannot fill SMILES |
| P69905-like heme b with Fe ligand part plus O2 | Two ligand groups; iron remains the heme ligand part, not a third ligand; evidence stays with the correct site |
| Generic substrate binding site and reaction text only | `no_matches` for identified ligands, but site and reaction sections are displayed |
| One ligand, failed CCD lookup | Core `success`, null structure, explicit `unavailable` warning; BindingDB unaffected |
| UniProt 404, BindingDB success | Both source rows saved; UniProt `Unavailable`; overall exit 1 |
| Both sources complete, no UniProt ligand annotations | UniProt count zero, no global absence claim, overall exit 0 |

No implementation work is complete until fixture-backed tests cover these cases, saved v4 data round-trips through the reader, all index links resolve to generated current-run files, and the complete existing suite passes after intentional schema/presentation updates.

## Tests and implementation sequence

Use checked-in minimal fixtures and HTTPX MockTransport; fixed UTC/monotonic clocks and no real sleeps or live tests in the normal suite.

- Entry identity: exact/secondary mappings, mismatches, reviewed checks, 404, empty-body retries, malformed entry, both sources failing independently, and aggregate exit codes.
- Extraction: multiple sites for one ligand, same ID across cofactors and binding features, different IDs with same name, name-only groups, generic substrates, unknown features, no annotations, and reactions with several chemical cross-references that must not become binding records.
- Coordinates: exact point/range, residue letters, uncertain endpoints, missing/unknown positions, alternate sequences, polyprotein containment/overlap, ligand parts, invalid coordinates, and evidence retained at annotation level.
- Structures: exact paired cross-reference validation, pagination, 204, wrong resource/chemical ID, stereo preference/fallback, ambiguous candidates, missing descriptors, enrichment outage, shared-ID cache, throttling/retries, and ATP's CHEBI:30616 versus CHEBI:15422 mismatch.
- Persistence: v1/v2/v3 BindingDB upgrades, v4 round trips and type/invariant validation, response manifests/checksums, no stale data after reruns, safe raw cleanup, write failures, and unchanged BindingDB raw handling.
- Presentation: index links resolve to generated files, two independent source statuses/counts, annotation-versus-measurement labels, conditional sections, full site/evidence retention, escaped remote text, selected SMILES provenance, no synthetic affinities, and TSV round trips.

Implement in this order: shared v4 models/reader; UniProt HTTP/extraction; optional CCD enrichment; source-specific renderers; independent orchestration and storage; README updates. Run relevant pytest tests after each code increment and the complete suite before completion. Preserve all current BindingDB behavior tests except presentation/schema assertions intentionally updated by this spec.

## Source references and verified contracts

- [UniProt sequence features](https://www.uniprot.org/help/sequence_annotation) and [UniProt manual](https://web.expasy.org/docs/userman.html): ligand/site annotation semantics and coordinates.
- [UniProt ligand annotation paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9825770/): ChEBI-based cognate ligands, chemical identity, and annotation evidence.
- [P00519 entry JSON](https://rest.uniprot.org/uniprotkb/P00519.json): observed binding-site, cofactor, and reaction JSON structures.
- [RCSB Search API](https://search.rcsb.org/) and [Data API](https://data.rcsb.org/): CCD identity lookup and descriptors.
- [CCD ATP JSON](https://data.rcsb.org/rest/v1/core/chemcomp/ATP): observed ChEBI mapping and stereo/non-stereo SMILES fields.

### GO annotation index section

The additive [GO annotation specification](go-annotations.md) defines UniProt GO term extraction, audited QuickGO relation enrichment, optional saved context, and MF/BP/CC presentation in `index.html`. This enrichment does not change ligand selection or counts.

### Accession-level function description in the index

UniProt supplies accession-level free-text function descriptions in `comments` entries with `commentType=FUNCTION` (see the [UniProt function guide](https://www.ebi.ac.uk/training/online/courses/uniprot-exploring-protein-sequence-and-functional-info/exploring-a-uniprotkb-entry/the-entry-view/function/)). The adapter already retains these in `source_context.context_comments`. Render their full, untruncated text before the GO section of `index.html`, with a UniProtKB/Swiss-Prot accession/function link and original retrieval date. Preserve any molecule/isoform label so qualified descriptions do not become general assertions. Escape descriptions and labels as HTML text. Exclude other comment types from this section. Missing FUNCTION text is reported as not supplied; a failed/absent UniProt source is reported as unavailable. Because UniProt offers this description, no `profunc` fallback is introduced. Rendering is offline and also applies to reused source results; no extra API request or schema change is required.

Function readability: retain separate source `texts[]` values and split embedded line endings into separate nonempty paragraphs in source order. For source paragraphs supplied as a single long block, add display-only paragraph boundaries at sentence-ending punctuation followed by whitespace and an uppercase letter, once at least 600 characters have accumulated. Avoid boundaries after common abbreviations (`e.g.`, `i.e.`, `Dr.`, `Mr.`, `Mrs.`, `vs.`, `Fig.`, `No.`). Keep all source wording and sentence order; disclose that additional paragraph breaks were added for readability. Do not truncate or summarize. Render explicit paragraph margins of 1.75em (a blank line between paragraphs), 1.75 line height, and a maximum reading width of 78ch. Keep attribution/retrieval metadata visually secondary, and retain isoform subheadings with clear spacing. Saved source data and provenance remain unchanged.
