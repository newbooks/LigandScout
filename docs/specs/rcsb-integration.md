# RCSB PDB ligand integration

Status: implemented. Add an independent RCSB PDB source to `ligandscout ACCESSION`, alongside BindingDB and UniProt. Keep the existing command and API-only architecture. This document describes experimental structure searches; existing UniProt-to-CCD enrichment remains separate.

## Scope and evidence

Find nonpolymer small molecules, cofactors, and ions with RCSB-reported contacts to a protein chain explicitly mapped to the query accession. Retrieve chemical identity and available SMILES from the Chemical Component Dictionary (CCD). Report contacting residues in PDB coordinates and map them to UniProt coordinates when an explicit alignment permits it.

Call these **structure-associated ligands**. A contact establishes proximity in a deposited structure, not measured affinity, physiological relevance, or binding to the isolated wild-type protein. Retain additives and ions with contacts, with their CCD names and any explicit bound flags retained in source neighbor rows; do not infer a biological role from a name. Exclude water components `HOH` and `DOD`, recording excluded counts. Do not equate missing data with absence of binding.

The first implementation covers experimental PDB entries and nonpolymer entities only. Exclude computed and integrative models, polymer partners, branched carbohydrates, and modified residues incorporated into a polymer. State these coverage limits in reports. Do not calculate new contacts, expand biological assemblies or crystal symmetry, infer covalent bonds from distance, or predict binding sites. Preserve neighbor `target_is_bound` values verbatim as source annotations, without interpreting undocumented values or requiring them for inclusion. Keep other optional role/covalent annotations in raw responses; no extra role-extraction algorithm is required in this increment.

## Implementation boundaries

The essential path is accession → exact protein entities → ligand instances with target contacts → CCD chemistry. UniProt residue mapping is an enrichment step: missing mapping must not hide a ligand with valid structural contacts. Keep source selection, coordinate mapping, and rendering as separate functions so each can be tested with small fixtures.

Implement this as one adapter plus source-specific schema/report helpers. Reuse the existing HTTP, storage, and record machinery. Do not introduce a graph client, coordinate parser, assembly builder, database mirror or concurrency layer. The shared `--force` option and completed-output reuse follow [query reuse](query-reuse.md). The schema migration is additive; the BindingDB and UniProt selection algorithms remain unchanged.

There is no request-count or ligand-count cap. Sequential requests can take minutes for well-studied proteins; log progress after each search page and entry, plus a concise completion summary. Do not print one terminal row per ligand. Keyboard interruption follows the existing CLI behavior (exit 130), rather than becoming a source failure.

## APIs and identity

No registration or API key is required. Use the existing synchronous HTTPX client and audited HTTP helper; no new client library or downloaded database is needed.

### 1. Search protein entities

Send `GET https://search.rcsb.org/rcsbsearch/v2/query` with a URL-encoded `json` parameter:

```json
{
  "query": {
    "type": "terminal",
    "service": "text",
    "parameters": {
      "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
      "operator": "exact_match",
      "value": "P00519"
    }
  },
  "return_type": "polymer_entity",
  "request_options": {
    "results_content_type": ["experimental"],
    "paginate": {"start": 0, "rows": 100},
    "sort": [{"sort_by": "rcsb_id", "direction": "asc"}]
  }
}
```

Use the normalized supplied accession, without sequence similarity or protein-name fallback. A successful UniProt primary/secondary accession resolution must supply a different primary accession as an additional exact search; preserve both queried IDs and their resolution evidence. Union and deduplicate entity IDs. If UniProt fails, still search RCSB with the original accession. Do not strip isoform suffixes or silently accept related accessions.

Advance `start` by the number of results returned until `total_count` is reached. Validate `result_type`, result identifiers, stable totals, page progress, and duplicate IDs within each search. Require at most 100 rows per page and exactly 100 for each nonfinal page; reject a page that exceeds the remaining total. A zero total requires an empty or absent result set; unexpected changes fail the source rather than claiming completeness. HTTP 204 on the first search page means no hits; on a later page it means an incomplete search. A valid JSON response with `total_count: 0` is also empty. No silent result cap. Sort identifiers locally for deterministic processing.

### 2. Retrieve and validate entities and entries

Use these Data API endpoints, caching successful objects by endpoint within one run:

| Object | Endpoint under `https://data.rcsb.org/rest/v1/core` | Purpose |
| --- | --- | --- |
| Polymer entity | `/polymer_entity/{entry}/{entity}` | Exact accession references, sequence alignment, target chain IDs, sequence and construct metadata |
| Entry | `/entry/{entry}` | Experimental method, title, resolution when supplied, citation, revision metadata and nonpolymer entity IDs |
| Nonpolymer entity | `/nonpolymer_entity/{entry}/{entity}` | CCD identity and instance IDs |
| Nonpolymer instance | `/nonpolymer_entity_instance/{entry}/{label_asym_id}` | Instance identity and `rcsb_target_neighbors` contacts |
| Chemical component | `/chemcomp/{ccd_id}` | Name, formula, SMILES and InChIKey |

Parse entity search identifiers with `rsplit("_", 1)` into nonempty entry and positive decimal entity IDs. Require alphanumeric entry IDs, preserving their complete length. Treat identifiers as URL-encoded path segments; do not assume every PDB ID will remain four characters. Check returned IDs and parent entry/entity relationships against the request. For polymer entities, require at least one reference object with `database_name == "UniProt"` and an accepted exact `database_accession` in `rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers`. Matching values in different array objects is insufficient. Reject mismatched entities into the selection audit. Require `rcsb_entry_info.structure_determination_methodology == "experimental"`. A different nonempty methodology excludes that entry with reason `unsupported_methodology`; missing or malformed methodology fails validation.

Collect all accepted target `asym_ids` (label chain IDs) and their entity IDs per entry. Enumerate `non_polymer_entity_ids` from entry container identifiers, then `rcsb_nonpolymer_entity_container_identifiers.asym_ids` on each nonpolymer entity. Missing optional collections mean empty as enumerated in the validation contract; malformed present collections are failures. Multiple target copies and multiple ligand copies remain distinguishable. Fetch an entry and each CCD component only once per run.

For example, the live `1IEP/1` polymer entity references **P00520**, not P00519. It must be excluded from a P00519 search unless a separate exact accepted reference exists. Never assign all ligands in a matching entry to every protein in that entry.

### 3. Select ligand contacts

Read `rcsb_target_neighbors` on each nonpolymer instance. Keep only neighbors whose `target_entity_id` and `target_asym_id` both identify an accepted target chain in that entry. Retain the source's `distance` in angstroms and contact fields including `atom_id`, `target_atom_id`, `target_comp_id`, `target_seq_id`, `target_auth_seq_id`, and `target_model_id` when present. Do not impose a second distance threshold or call the source's precomputed neighborhood a newly calculated binding pocket.

A ligand instance is included only if it has at least one validated contact to an accepted target chain. Nonwater instances with nonempty neighbors but no accepted target contact are excluded with reason `no_target_contact`; absent or empty neighbors use `no_reported_contacts`; water is excluded with reason `water`. An absent contact collection means no reported contacts, not proof of no physical contact. Report that limitation and count these instances separately. A malformed contact collection or a failed required instance request fails the source; never convert an unavailable contact service to a successful zero-match result.

Group contacts by ligand instance, target chain, and target model. Deduplicate residues by target label sequence position within that group, preserving all returned source neighbor rows separately (these are not an exhaustive list of atom pairs). A ligand contacting several target chains retains separate sites. Do not synthesize alternate locations, insertion codes, ligand model IDs, assembly IDs, or author-to-label mappings that the response does not provide. Unavailable fields remain null.

## Binding-site coordinates

Use `rcsb_polymer_entity_align` entries with `reference_database_name == "UniProt"` and an accepted exact accession. Each `aligned_regions` segment maps entity position `p` using:

```text
entity_beg_seq_id <= p < entity_beg_seq_id + length
uniprot_position = ref_beg_seq_id + p - entity_beg_seq_id
```

These are 1-based coordinates. Validate positive segment starts/lengths and entity bounds. Do not extrapolate beyond segments. Map only when all covering segments for that accession agree; conflicting or absent mappings leave the UniProt position null with an explicit reason. An absent alignment does not discard a valid structural contact. Malformed alignment data produces a warning and unavailable mapping, preserving raw data and PDB sites.

Each site retains entry ID, target entity, label chain, target model, residue names, label sequence positions, author residue positions as supplied, contact distances and mapping provenance. Target author chain remains null in this increment, as specified below. Never use author residue numbers directly as UniProt positions. Preserve insertion-code information if supplied; do not manufacture it from numeric author positions.

If the resolved UniProt sequence is available, check mapped positions against its bounds and show its residue alongside the observed PDB residue. Record substitutions/mismatches and available mutation/construct metadata; do not assume an accession-mapped construct is wild type. An out-of-bounds position leaves that residue unmapped with a warning. A residue substitution alone does not invalidate a valid positional alignment: retain the mapped position and label the sequence mismatch. If UniProt is unavailable, retain RCSB mappings labeled as not independently sequence-checked. Chimeras and processed proteins require segment-specific mappings, not a global offset.

Report sites as residue bullet lists, grouped by entry/chain/model, with PDB coordinates always available when supplied. Display UniProt coordinates separately and explain unmapped residues. Say “RCSB-reported residue contacts”; do not imply that every intervening residue in a range contacts the ligand. Contact scope is the source's reported neighborhood; assembly/symmetry contacts and all conformers/models are not guaranteed covered.

## Chemical identity and SMILES

Use each included nonpolymer instance's CCD ID, cross-checked against its parent entity. Fetch that exact component directly; no ChEBI search is necessary. Reuse the CCD request/descriptor parsing code where practical, but separate direct CCD retrieval from UniProt's exact ChEBI-to-CCD selection policy.

Require matching `rcsb_id` and `chem_comp.id`. Prefer nonempty `rcsb_chem_comp_descriptor.SMILES_stereo`; otherwise use `SMILES` and label it nonstereospecific. Preserve InChIKey when available. Record provider, CCD ID, chosen descriptor field, request ID and retrieval time. Do not generate, neutralize, or merge structures, or substitute a similarly named component. Missing descriptors, ions without descriptors, and failed CCD enrichment retain the ligand with a stated reason and null SMILES. CCD failure is an enrichment warning rather than failure of the structural search.

Group the same exact CCD ID into one ligand record across entries, retaining every distinct instance and site. Different CCD IDs remain separate even if names or SMILES match. Do not merge RCSB records with BindingDB or UniProt records; cross-source agreement is not independent confirmation of affinity.

## Saved-result contract

Add source `RCSB PDB`, slug `rcsb`, policy `rcsb-structural-contacts-v1`. Extend the normalized schema to **version 5** for structural records and counts, updating readers and writers together. Continue reading BindingDB v1–v4 and UniProt v4 without changing their meaning. Keep the complete existing TSV column set and order unchanged. Use the existing `binding_sites`, `structure_provenance`, and `source_values` JSON-valued columns for structural details; add no columns in this increment.

Reuse `Measurement` as the shared record container:

- `compound_id`: `CCD:{id}`; name from `chem_comp.name`, falling back to `pdbx_entity_nonpoly.name`, then CCD ID. For conflicting entity names without a CCD name, choose the first by sorted entry/entity ID and retain original entity names in `source_context.ligand_entities`, keyed by `entry_entity`.
- `evidence_type`: `structure_associated_ligand`; all affinity/assay fields and top-level PMID/DOI remain null. Do not import `rcsb_binding_affinity`, which may duplicate BindingDB evidence. Never create ligand records from `rcsb_target_cofactors`, `rcsb_related_target_references`, chemical similarity, or the presence of a ligand anywhere else in the entry: only the validated instance-neighbor path establishes inclusion. Entry citations are retained in the entry context referenced by observations, not as fabricated measurements.
- `source_url`: `https://www.rcsb.org/ligand/{id}`; `target_accessions`: sorted exact accepted accession references actually present on included target entities; `target_match_status`: `exact_uniprot_reference`.
- `binding_sites`: structural site objects containing the coordinates and provenance defined above; keep UniProt annotation site validation source-specific.
- `structure_provenance`: direct CCD lookup outcome and descriptor provenance.
- `source_values.observations`: entry, ligand entity/instance, target chain/model groups, raw contact fields and request IDs. Entry and entity metadata are stored once in source context and joined through these IDs. Stable observation key is entry + ligand label chain + target label chain + target model (null explicitly represented).

Use `counts.record_kind = "structural_ligand"`. Add nullable structural counters: `candidate_entries`, `accepted_target_entities`, `examined_ligand_instances`, `included_ligand_instances`, `excluded_water_instances`, `instances_without_reported_contacts`, `instances_without_target_contacts`, `contact_sites`, `contact_residues`, and `mapped_contact_residues`. Define `contact_residues` as unique residue positions per site, and `contact_sites` as ligand-instance/target-chain/model groups. `unique_ligands` is the number of CCD groups; SMILES counters count those groups. Measurement and UniProt-only annotation counters are null for RCSB; structural counters are null for other sources. The exact definitions and invariants below are required, not deferred design choices.

`selection_audit` retains accepted/rejected entities and examined/excluded instances with reasons and request IDs. `source_context` retains queried accessions, resolution evidence, coverage limits, and entry/entity metadata. Missing-contact totals belong in `counts`; mapping warnings belong in `warnings`. Unknown entry revision/version fields remain null; record retrieval times regardless.

Save alongside existing outputs:

```text
results/ACCESSION/
  sources/rcsb.json
  sources/rcsb-requests.json
  sources/rcsb-raw/000001.json
  rcsb-ligands.html
```

Reuse the UniProt request-manifest format, body hashes, numbered raw response files and non-JSON body handling. Persist every attempt, redirect and response, including failed CCD calls. Extend staged overwrite and obsolete-owned-raw cleanup to RCSB; do not delete other sources' files. Reruns replace stale RCSB rows after empty or failed searches. No directory-wide atomicity or concurrent same-accession run guarantee is introduced.

## Orchestration, failures and presentation

Run BindingDB, UniProt and RCSB independently in the existing workflow. Reuse successful UniProt resolution/sequence as optional RCSB context, without requiring UniProt success. Sequential requests use the shared audited helper: at least 0.5 seconds between request starts, 30-second operation timeouts, three attempts, existing backoff/Retry-After limits and HTTPS redirect rules. Add explicit support for first-page RCSB search HTTP 204; do not treat Data API 404 or empty HTTP 200 bodies as empty results. Preserve all queries and page offsets. Use separate source-local audits and caches. Do not share CCD responses across UniProt and RCSB in this increment; this keeps each manifest self-contained. Within RCSB, cache each successful endpoint response and link every use to its original request ID.

A complete search with included ligands is `success` (with `error: null`); a complete search without included ligands is `no_matches`. Explain whether there were no exact structures, no eligible entities, or no reported nonwater target contacts. Search/pagination, entry, target identity, entity or contact retrieval failures make RCSB `failed` and the command exit 1, while other sources still finish. On a required failure, stop RCSB processing, return `records: []`, set `counts.complete` false and all numeric counters null. Keep completed selection audit/context and all raw exchanges for diagnosis; do not publish partial ligand rows. If failure occurs during CCD enrichment, it is optional and follows the enrichment rule instead; the core-failure rule applies only to required calls and validation. This preserves the existing status/record invariant. Optional descriptor and mapping failures retain structural evidence with warnings. Do not silently skip unavailable structures.

Add RCSB to `index.html`, `report.md`, and the combined TSV. Audit all source dispatch branches in `report.py`, including `source_slug`, `compound_url`, `evidence_note`, `render_summary`, HTML and Markdown: a generic “not UniProt” branch must not treat RCSB as BindingDB. The compact RCSB summary reports unique structural ligands and contact sites; failed counts read “Unavailable”. Keep the index Records column as grouped ligand-record count and label its RCSB meaning in the source page. Link `rcsb-ligands.html` from the index; keep the terminal compact and do not list individual ligand-page filenames. Label RCSB records as structural ligands/sites, not measurements. Do not sum unique ligands across sources.

The source HTML table has ligand/CCD link, available SMILES, supporting structures, and contact-site bullet lists. Link structures to `https://www.rcsb.org/structure/{entry}` and provide Data API provenance links. Show experimental method/resolution only when supplied, target chain identity, mapping warnings and coverage limits. Group long per-structure evidence in expandable details. Use plain short text, escaped remote content, selectable SMILES, and links to saved JSON; do not render JSON objects in HTML. Keep source evidence distinct from physiological interpretation.

## Concrete adapter and validation contract

Implement `fetch_rcsb(accession, *, client, sleep, now, monotonic, uniprot_context=None) -> SearchOutcome`. It normalizes the accession and creates its own `Audit`. Set `raw_body=None`; all RCSB bodies are in `exchanges`. `run_search` calls it after UniProt, then returns sources in order BindingDB, UniProt, RCSB PDB. Pass `uniprot.result.source_context["entry"]` only when UniProt status is `success` or `no_matches` and that entry is present. Do not read saved files or issue another UniProt request. The existing context path is `source_context["entry"]`; cover its absence and the zero-ligand UniProt case in tests.

Accepted accessions are the original query plus the resolved primary only when the validated entry has that query as primary or secondary. Do not add all other secondary accessions. Set `preferred_accession` to the resolved primary, otherwise the query. Search original first, primary second if different. RCSB target references must match one of this set; an unrelated reference is an audited exclusion, not an identity-validation error.

Process all candidate polymer entities before ligand enumeration so every accepted target chain in an entry is known. Fetch entries only for entries containing an accepted entity. Once an entry is accepted, fetch every listed nonpolymer entity and instance in sorted order. Entry `non_polymer_entity_ids` absent means no enumerated nonpolymer entities; its presence requires a list of unique positive decimal ID strings. Water entities may be absent from this list: water counters describe enumerated instances only, not all crystallographic waters. Do not fetch water entities through unrelated `entity_ids`.

### Required fields and missing data

Validate only fields consumed by this adapter; preserve unknown fields in raw responses. Object/list mismatches and invalid required types are core schema failures. Reject booleans where integers/numbers are required, and reject nonfinite distances. Empty optional strings become null. Lists of IDs must contain unique nonempty strings. Require:

| Object | Required identity/data |
| --- | --- |
| Search page | `result_type == "polymer_entity"`, nonnegative integer `total_count`, list `result_set` (may be omitted only for zero total); each row has a valid `identifier` |
| Polymer entity | matching `rcsb_id` and container `entry_id`/`entity_id`; nonempty `asym_ids`; `reference_sequence_identifiers` list (absent means empty); positive `entity_poly.rcsb_sample_sequence_length` |
| Entry | matching `rcsb_id` and container `entry_id`; methodology as above |
| Nonpolymer entity | matching `rcsb_id` and container parents; nonempty `asym_ids`; nonempty `nonpolymer_comp_id`; `pdbx_entity_nonpoly.comp_id`, when present, must agree |
| Nonpolymer instance | matching `rcsb_id` (`entry.label_chain`), container parents, `asym_id`, `comp_id` matching entity CCD ID |
| Neighbor | nonempty `target_entity_id` and `target_asym_id`; for accepted targets, positive integer `target_seq_id`, nonempty `target_comp_id`, finite nonnegative numeric `distance` |

For neighbor rows on other chains, require identity fields to permit exclusion but do not require their residue fields. A label chain belonging to a known accepted entity but paired with a different entity ID is a schema failure. A chain absent from the accepted chain map is an unrelated target and is excluded, even if its entity ID matches an accepted entity; do not infer chain membership from entity ID alone. Duplicate label chain IDs assigned to different accepted entities fail validation. For accepted neighbors, optional `target_model_id` must be a positive integer, `target_auth_seq_id` an integer (zero/negative author numbering is allowed), and supplied atom/alternate/bound flags strings. Missing atom names do not invalidate a residue contact. A present neighbor `comp_id` must agree with the ligand. Require target label sequence positions within the target entity length. Missing/empty `rcsb_target_neighbors` gives exclusion `no_reported_contacts`; a nonempty collection containing only other targets gives `no_target_contact`. Null collections are malformed, not empty.

Do not zip polymer entity `asym_ids` and `auth_asym_ids`: correspondence is not guaranteed by list position. Target author chain is null in this increment; the label chain plus source author residue number is sufficient and clearly labeled. Ligand author chain/position come directly from the instance container (`auth_asym_id`, `auth_seq_id`, both optional strings). No polymer-instance request is needed. Polymer name/construct context comes from `rcsb_polymer_entity.pdbx_description`, `rcsb_polymer_entity.pdbx_mutation` when supplied, and `entity_poly.pdbx_seq_one_letter_code_can`; retain canonical sequence text as source metadata without using its formatted length for indexing.

Optional entry fields are `struct.title`, `exptl[].method`, `rcsb_entry_info.resolution_combined` (list of positive finite angstrom values), `rcsb_primary_citation` (title/DOI/PubMed ID), and `rcsb_accession_info` (revision date, major/minor revision). Missing fields remain null/empty; malformed optional metadata emits `invalid_optional_metadata` and is omitted. `source_version` is null because multiple entries have distinct revisions; save revisions per entry. `retrieved_at` is the last completed RCSB request timestamp (including enrichment). `request` contains `search_endpoint`, ordered `query_accessions`, `page_size: 100`, `return_type: "polymer_entity"`, and `results_content_type: ["experimental"]`.

### Mapping decisions

Collect alignments only for exact accepted references. Choose the preferred accession if any of its valid segments covers a residue; otherwise choose the original query if it has coverage. If preferred-accession covering segments conflict, return `conflicting_alignment`; do not fall back to the query to conceal that conflict. Do not combine positions from different accessions into one coordinate system. Multiple covering segments for the chosen accession must agree. A malformed matching alignment invalidates all mappings for that entity with `invalid_alignment`; unrelated malformed optional alignments are ignored. Segment bounds use `entity_poly.rcsb_sample_sequence_length`.

Reason precedence is: invalid matching alignment → `invalid_alignment`; no matching alignment objects → `no_alignment`; matching objects but no segment covering the residue → `outside_alignment`; inconsistent covering positions → `conflicting_alignment`; mapped primary position outside the available sequence → `uniprot_out_of_bounds`; otherwise `mapped`. Missing or empty `aligned_regions` on a matching alignment is invalid, whereas an absent/empty root alignment collection means no alignment. Null or nonlist root alignment data invalidates mapping for the entity. Mapping reasons are `mapped`, `no_alignment`, `outside_alignment`, `conflicting_alignment`, `invalid_alignment`, or `uniprot_out_of_bounds`. Sequence check is `match`, `mismatch`, `unknown_residue`, or `not_checked`. Check only when the mapped accession equals the available sequence's primary accession. Convert the standard 20 three-letter residue codes to one-letter codes with a fixed lookup; other components are `unknown_residue`. A mismatch retains the position and emits a warning. No residue realignment or mutation correction is performed.

### Exact structural objects

Each `binding_sites` object has these keys (null is explicit for unavailable scalar fields):

```text
entry_id: str
ligand_entity_id: str
ligand_asym_id: str
ligand_auth_asym_id: str|null
ligand_auth_seq_id: str|null
target_entity_id: str
target_asym_id: str
target_auth_asym_id: null
target_model_id: int|null
residues: list[Residue]
request_ids: list[int]
```

`Residue` has `label_seq_id: int`, `auth_seq_id: int|null`, `comp_id: str`, `distance: number` (minimum returned distance for this residue), `uniprot_accession: str|null`, `uniprot_position: int|null`, `mapping_reason: str`, `uniprot_residue: str|null`, `sequence_check: str`, and `alignment_segments: list[dict]` (the source segments supporting the mapping). Set accession and position only for a valid mapping; otherwise both are null and segments empty. A valid mapping to a secondary accession is labeled with that accession and never relabeled as primary. `uniprot_residue` is null unless the sequence check can read that primary sequence position. Preserve insertion/alternate information only in the original neighbor rows; no insertion-code field is inferred. For repeated rows at one label position within a site, residue names must agree. Merge missing author positions with a single supplied value; two different nonnull author positions fail validation. Do not let optional-field absence alone invalidate otherwise consistent contacts.

Each `source_values.observations` object has `entry_id`, `ligand_entity_id`, `ligand_asym_id`, `target_entity_id`, `target_asym_id`, `target_model_id`, `neighbor_rows` (all matching rows verbatim, including duplicate rows), and `request_ids`. Store entry metadata in `source_context.entries` keyed by entry ID, target construct context in `.target_entities` keyed by `entry_entity`, and nonpolymer names/CCD IDs in `.ligand_entities` keyed by `entry_entity`. There is exactly one observation and one site per unique `(entry_id, ligand_asym_id, target_asym_id, target_model_id)` tuple. Sort records by CCD ID, sites/observations lexically by the first three tuple fields then model (null before integers), and residues by label position. Keep source row order within `neighbor_rows`.

Record `source_record_id` and `target_id` are null because each grouped ligand can cover many entries/entities. `target_description` is null; per-entity names remain in context. `ligand_roles` and `annotation_evidence` are empty lists; do not manufacture annotation evidence or classifications. All other unpopulated `Measurement` fields retain dataclass defaults.

Direct CCD `structure_provenance` has exactly `provider: "RCSB CCD"`, `status` (`available`, `no_descriptor`, `unavailable`), `ccd_id`, `smiles_field` (chosen descriptor key or null), and `request_ids`. Missing SMILES with a valid InChIKey retains the key. Failed/mismatched/malformed components set both descriptors null and `unavailable`; warn once per CCD ID. Do not reuse UniProt's ChEBI provenance shape or its annotation validator.

### Counts, audits and schema migration

For complete RCSB results define:

| Counter | Definition |
| --- | --- |
| `candidate_entries` | Unique entry IDs in the union of entity search results, before reference validation |
| `accepted_target_entities` | Exact-reference entities in experimentally validated entries |
| `examined_ligand_instances` | Enumerated nonpolymer instances whose required response was validated |
| `included_ligand_instances` | Unique `(entry, ligand label chain)` pairs having at least one accepted site |
| `excluded_water_instances` | Enumerated instances excluded by CCD ID HOH/DOD |
| `instances_without_reported_contacts` | Nonwater instances with absent/empty neighbor collection |
| `instances_without_target_contacts` | Nonwater instances with neighbors but no accepted target contact |
| `contact_sites` | Number of site objects across all records |
| `contact_residues` | Sum of residue-list lengths across sites |
| `mapped_contact_residues` | Residues with a nonnull UniProt position, counted separately per site |
| `unique_ligands` | Number of grouped CCD records |
| `smiles_available`, `smiles_unavailable` | CCD records with and without nonnull SMILES |

Require `examined_ligand_instances == included_ligand_instances + excluded_water_instances + instances_without_reported_contacts + instances_without_target_contacts`; `mapped_contact_residues <= contact_residues`; `unique_ligands <= included_ligand_instances`; and SMILES counts sum to `unique_ligands`. Water exclusion takes precedence over contact checks, after instance identity validation. Zero results have numeric zeros for applicable counters, not nulls. All counters use strict nonnegative integers.

The source-context root has exactly `query_accessions`, `preferred_accession`, `uniprot_resolution`, `coverage_limits`, `entries`, `target_entities`, and `ligand_entities`. `query_accessions` preserves search order; `uniprot_resolution` is the validated query/primary pair or null. `coverage_limits` is a list of the scope limitations above. Entry objects retain the optional entry fields under their original API names; target entity objects retain exact reference objects, entity length, optional construct fields and alignment objects; ligand entity objects retain CCD ID and available original name. Mapping warnings stay in `warnings`, and missing-contact counts stay in `counts`, avoiding duplicate totals. Keep all root keys present on failure, with empty dictionaries for objects not yet retrieved. Metadata objects have the specified selected source fields plus `request_ids`; source-context dictionary keys and observation IDs must agree. Validate these relationships and count invariants in `rcsb_schema.py`.

For site/observation `request_ids`, record the sorted unique logical request IDs for the validating entry, polymer entity, nonpolymer entity and instance responses. Component request IDs belong in structure provenance. Cache hits reuse those same IDs. Source-context objects also carry their own `request_ids`; IDs refer to the local manifest and must exist in it.

Schema-v5 counts have the union of the existing v4 keys and these ten structural counters. Upgrade old results to v5 in memory with structural counters null; their existing values remain unchanged. Native BindingDB/UniProt adapters may still produce older schemas internally; shared upgrade runs before serialization. For BindingDB and UniProt validation, create a v4 projection (`schema_version=4`, v5-only counters removed) and run their existing validation rules; do not call the public reader recursively. Restore the v5 result afterward. The upgrade function must be idempotent: upgrading a v5 object returns it unchanged and does not enter the BindingDB legacy-upgrade path. RCSB is supported only at v5 and uses its own validator (`rcsb_schema.py`). For RCSB the activity counters and five annotation-only counters are null; the shared SMILES counters are numeric. `status == "success"` iff records are nonempty, `no_matches` requires empty records and no error, and `failed` requires empty records, an error and `complete: false`. Never permit failed-source records for any source. New readers accept native v5 for all three sources; legacy reader tests remain required.

RCSB `selection_audit` has exactly `outcome_reason`, `targets`, and `instances`. Target rows have `entry_id`, `entity_id`, `included: bool`, `reason: str|null`, `request_ids: list[int]`; instance rows have `entry_id`, `entity_id`, `asym_id`, `ccd_id`, `included: bool`, `reason: str|null`, `request_ids: list[int]`. Included rows use `reason: null`. Allowed exclusion reasons are `accession_mismatch`, `unsupported_methodology`, `water`, `no_reported_contacts`, `no_target_contact`. On required failure set outcome `incomplete`; otherwise `no_exact_structures` for zero search hits, `no_eligible_entities` for no accepted entities, `no_reported_target_ligands` for accepted entities without included ligands, and null for success. Do not add fabricated audit rows for requests that never completed validation.

Warnings retain the existing `{code, message, role}` shape. Use request roles `rcsb_search`, `rcsb_entry`, `rcsb_polymer_entity`, `rcsb_nonpolymer_entity`, `rcsb_nonpolymer_instance`, and `rcsb_component`. `chemical_id` is the CCD ID for component calls and null otherwise. HTTP helper allows 204 for `rcsb_search`; the adapter enforces page position. Schema/identity errors use `SourceError(kind="schema", message=..., http_status=...)`; transport/HTTP errors preserve the helper's error. Catch expected validation/request errors inside RCSB, not arbitrary programming exceptions.

No cross-source shared cache or global request numbering is introduced. Each RCSB logical call has one local request ID across retries/redirects; each exchange gets a unique sequence. Persist the complete source-local manifest even on failure. Extract reusable manifest writing in `storage.py` parameterized by source slug; preserve UniProt filenames/content and BindingDB raw behavior. Local report links resolve under the accession output directory. Manifest `raw_path` values resolve relative to `sources/`, following the existing UniProt convention. External provenance links remain absolute HTTPS URLs.

## Implementation sequence

1. Add schema-v5 upgrade/reader support and the RCSB validator, retaining legacy fixtures.
2. Implement accession/entity selection and instance contact extraction with source-local audit; test the baseline fixture before adding optional enrichments.
3. Add direct CCD descriptors and coordinate mapping as separate enrichment functions. CCD failures never change structural selection.
4. Add orchestration, manifest storage and source-specific report dispatch. Verify overwrite and failed-result behavior.
5. Run the relevant and full pytest suites, then update the README and existing specs for the third source and schema version.

## Implementation and acceptance tests

Add `rcsb.py` for search/extraction/mapping and `rcsb_report.py` for source HTML. Extend models/schema validation, saved-result readers, storage, source dispatch, orchestration, CLI summary, fixtures and relevant specifications. Extract shared CCD retrieval/descriptor helpers without changing UniProt chemical-match rules.

Use a small, hand-verifiable end-to-end fixture as the baseline: one accepted target entity with chains A/B; four enumerated nonpolymer instances (one STI contacting A and B, one CL contacting A, one water, one molecule contacting only an unrelated chain). This yields 1 candidate entry, 1 accepted entity, 4 examined instances, 2 included instances, 1 excluded water, 0 without reported contacts, 1 without target contacts, 2 unique ligands and 3 sites. Give each site one residue and leave one residue outside alignment: 3 contact residues and 2 mapped residues. Supply SMILES only for STI: available/unavailable counts are 1/1. Assert the saved JSON, TSV, HTML and reader agree with these values. Add a variant with absent neighbors to test the separate no-reported-contact counter.

Use mocked pytest fixtures; routine tests must not depend on live services. Cover:

- Exact accession search, primary/secondary resolution and independent operation after UniProt failure; exclude a P00520-only entity for a P00519 query.
- Multiple pages, zero hits/first-page 204, later-page 204, changing totals, duplicate IDs, empty/malformed bodies, retries, redirects, throttling and failed required requests.
- Multiple proteins per entry: only contacts to verified target entity/chain pairs are included; test multiple copies, waters, ions, no contacts, unrelated-chain contacts and malformed contacts.
- Multiple structures sharing a CCD ID and distinct ligand copies/sites; deterministic grouping, counts, ordering, request reuse and preserved provenance.
- Stereo/fallback/missing SMILES, CCD identity mismatch and enrichment outage without loss of structural records.
- Multiple alignment segments, boundaries, gaps, conflicting mappings, chimeras, author/label differences, missing mappings/model fields, sequence mismatches and UniProt unavailability.
- Schema-v5 round trips and legacy readers; mixed-source TSV, source-specific counts, failure totals, escaped HTML, site bullets, working index links and no inline JSON.
- Overwrite after success, zero matches and failure, raw manifests/hashes, obsolete RCSB raw cleanup and preservation of unrelated source artifacts.

Run affected pytest tests and the full suite before completion. An optional manual API smoke check should verify one mapped protein with nonwater contacts and one zero-hit case; record observations without hard-coding mutable database counts. Update README and mark this spec implemented only after integration and checks pass.

## API references

- [RCSB Search API](https://search.rcsb.org/) and [search attributes](https://search.rcsb.org/structure-search-attributes.html).
- [Data API: endpoints, hierarchy, schemas and error handling](https://data.rcsb.org/index.html).
- Live response examples checked while drafting: [1IEP polymer entity 1](https://data.rcsb.org/rest/v1/core/polymer_entity/1IEP/1), including explicit UniProt alignment; [1IEP nonpolymer instance C](https://data.rcsb.org/rest/v1/core/nonpolymer_entity_instance/1IEP/C), including target-neighbor atom distances and residue identifiers. These illustrate API fields, not a P00519 binding claim.

Specification validation: checked the Search API query (including `rcsb_id` sorting), Data API schema fields, and live 1IEP nonpolymer entity/instance responses on 2026-11-05. The search response is mutable; its result count is not an acceptance-test constant. See [nonpolymer instance schema](https://data.rcsb.org/rest/v1/schema/nonpolymer_entity_instance) and [polymer entity schema](https://data.rcsb.org/rest/v1/schema/polymer_entity).

## Implementation verification

The adapter is implemented in `rcsb.py`, with structural validation in `rcsb_schema.py` and source-specific HTML/Markdown in `rcsb_report.py`. The shared CCD descriptor selector is reused without changing UniProt chemical matching. Endpoint responses are cached within the RCSB adapter; completed saved source results can also be reused across invocations under the shared [query reuse policy](query-reuse.md). Schema v5 is written for all three sources; legacy source contracts are projected and validated before upgrade.

Live API smoke checks on 2026-11-05: P00921 completed with 15 candidate entries, 14 CCD ligand groups, 113 contact sites, and SMILES for all 14 groups; P00699 completed with zero hits. These counts are observations, not test constants. Source reports were saved under `/tmp/ligandscout-rcsb-smoke/`. The positive run exercised transient retries and retained 451 mapped contact residues; coordinates were not independently checked against a UniProt sequence in this adapter-only run.

### CCD identifier links

In `rcsb-ligands.html`, each displayed `CCD:ID` links to `https://www.rcsb.org/ligand/ID`. URL-encode the component ID as a path segment and HTML-escape the displayed identifier. Keep the existing compound-name and CCD API provenance links.

### Residue mismatch presentation

Retain mapped residue mismatches in saved JSON and detailed contacts. Suppress their redundant `residue_mapping` warnings (`mapped, mismatch`) from terminal, index, Markdown and source-page summaries, without suppressing other mapping warnings. Highlight each mismatching contact in the HTML table in amber, retaining the explicit `sequence mismatch` label and explaining the color above the table. Matching contacts retain normal styling.
