# ChEMBL integration specification

Status: implemented. API contracts checked November 4, 2026. The CLI now searches BindingDB and ChEMBL and writes the multi-source outputs defined here.

This specifies a second simple-mode source for `ligandscout`, using only online
ChEMBL APIs. It extends [CLI and BindingDB](cli-bindingdb.md), retaining that
adapter's behavior. Where multi-source interfaces below differ from the existing
single-source contract, this specification takes precedence for the integration
increment. BindingDB request, parsing, retry, and inclusion policies remain
unchanged. Scope is occasional batches of roughly 20 accessions, invoked individually;
no batch command, database dumps, local mirrors, or similarity searches.

## Research and verified API behavior

Official references:

- [Data web services](https://chembl.gitbook.io/chembl-interface-documentation/web-services/chembl-data-web-services): resources, JSON, filters, pagination, and accession lookup.
- [Interactive API specification](https://www.ebi.ac.uk/chembl/api/data/docs).
- [Target questions](https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/target-questions): target relationships, variants, isoforms.
- [Assay/activity questions](https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/chembl-data-questions): assay types, confidence, standardization, and quality flags.
- [Interface URLs](https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/chembl-interface-questions): current report-card routes.

Read-only probes of the [release endpoint](https://www.ebi.ac.uk/chembl/api/data/status.json),
[exact P0DTC1 lookup](https://www.ebi.ac.uk/chembl/api/data/target.json?target_components__accession=P0DTC1),
and [example assay](https://www.ebi.ac.uk/chembl/api/data/assay/CHEMBL674637.json)
on November 4, 2026 returned the following observations. These
are smoke-check evidence, not permanent assertions about database coverage:

| Probe | Observation |
| --- | --- |
| `status.json` | `status: UP`, `chembl_db_version: ChEMBL_37`, release date `2026-05-01` |
| `target.json?target_components__accession=P0DTC1` | `targets: []`, `page_meta.total_count: 0` |
| Same lookup for `P0DTD1` | Includes `CHEMBL4523582`, Replicase polyprotein 1ab, and protein-interaction targets |
| Same lookup for `P00533` | Includes `CHEMBL203`, SINGLE PROTEIN, plus families, a chimera, and protein-interaction targets |
| `activity.json?target_chembl_id=CHEMBL203&order_by=activity_id&limit=1` | Activity `32260`, assay `CHEMBL674637`, molecule `CHEMBL68920`, document `CHEMBL1134862` |
| `assay/CHEMBL674637.json` | Confidence `8`, relationship `H`, despite the activity's exact target ID and binding assay type `B` |
| Batched `assay`, `molecule`, `document` collections | Each accepted its resource ID `__in` filter and `order_by` on that ID; two-ID assay and molecule batches also returned both requested entities; JSON collection keys were `assays`, `molecules`, `documents` |
| `assay/CHEMBL676150` and its activities | Direct assignment `D`, confidence `9`, type `B`; activity `1102403`, molecule `CHEMBL29197`, standardized IC50 `0.029 nM` |
| `assay.json?target_chembl_id=CHEMBL203&variant_sequence__isnull=false&limit=1` | Assay `CHEMBL939338` exposes `variant_sequence.accession=P00533`, `mutation=L858R`, `isoform=1` |

The homolog-assigned assay observation is a regression case: accession lookup alone does not exclude
assays mapped from homologs. A target name or activity target ID alone is insufficient.
P0DTC1 and P0DTD1 must not be substituted for each other, even when a processed
protein or domain overlaps. An empty P0DTC1 lookup means no exact coverage under
this policy; it does not establish that related viral proteins have no ligands.

## Command and inclusion policy

Keep the command syntax, normalization, validation, logging options, help behavior,
and explore placeholder. After this increment, default/explicit simple mode queries
BindingDB and ChEMBL, sequentially and independently. No new CLI options are needed
for this increment. Continue the other source after an expected service failure.

ChEMBL policy `exact-single-protein-v1`:

1. Resolve targets by exact component accession, not name or function. Inspect all
   returned targets; accept only `target_type == "SINGLE PROTEIN"` with exactly
   one component, with `component_type == "PROTEIN"`, whose accession equals the
   normalized query (case-sensitive equality after query normalization). Retain
   every eligible target ID, rather than arbitrarily selecting the first.
2. Inspect each assay associated with candidate activities. Include only assays
   whose target ID is eligible, `relationship_type == "D"`, `confidence_score == 9`,
   and `assay_type` is `B` or `F` (binding or functional). Require activity and assay
   target IDs and any nonnull assay types to agree; contradictory records fail
   validation. Assay metadata is authoritative when activity metadata is absent.
3. Exclude families, complexes, interactions, chimeras, homolog assignments, lower
   confidence assignments, and other assay types from the counted ligand set.
   Missing/null confidence or relationship is unverified and excluded, not guessed.
4. Retain variants under their mapped parent accession and label their supplied
   mutation/accession metadata. Null variant metadata means "not annotated", not
   proof of wild-type sequence. Domains, constructs, and isoforms may share a
   target accession; preserve descriptions and do not claim full-length identity.
5. Include all activity endpoints in eligible assays, quantitative or qualitative,
   without potency, activity-type, molecular-weight, molecule-type, approval, or
   positive/negative activity filters. Keep inactive and inconclusive records and
   quality/duplicate flags visible. A count means unique **tested compounds with
   eligible activity records**, not proven binders or a list of potent inhibitors.

This is a conservative target-assignment policy, not independent sequence or
experimental verification. ChEMBL's `D` denotes direct target assignment; it does
not by itself prove physical ligand binding. Do not apply these ChEMBL-specific
rules to BindingDB's descriptive labels. Its target identity remains unverified.
Future homolog, complex, wild-type-only, or potency modes require separate policy
specifications and explicit report labels.

## Online request pipeline

Base URL: `https://www.ebi.ac.uk/chembl/api/data/`. Use the existing injectable
HTTPX client, TLS verification, JSON Accept header, and versioned User-Agent.
No API key or new runtime client dependency is required. Base URL is a code
constant, not a user-supplied endpoint in this increment. Pass filters through
`params`; do not concatenate unescaped query strings.

| Stage | Relative endpoint and initial parameters | Envelope |
| --- | --- | --- |
| Release | `status.json` | Object containing `chembl_db_version` |
| Resolve accession | `target.json`, `target_components__accession=ACCESSION`, `limit=100`, `offset=0`, `order_by=target_chembl_id` | `targets`, `page_meta` |
| Activities per eligible target | `activity.json`, `target_chembl_id=CHEMBL_ID`, `limit=100`, `offset=0`, `order_by=activity_id` | `activities`, `page_meta` |
| Assay details | `assay.json`, `assay_chembl_id__in=ID1,ID2,...`, `limit=100`, `offset=0`, `order_by=assay_chembl_id` | `assays`, `page_meta` |
| Molecule enrichment | `molecule.json`, `molecule_chembl_id__in=ID1,ID2,...`, `limit=100`, `offset=0`, `order_by=molecule_chembl_id` | `molecules`, `page_meta` |
| Document enrichment | `document.json`, `document_chembl_id__in=ID1,ID2,...`, `limit=100`, `offset=0`, `order_by=document_chembl_id` | `documents`, `page_meta` |

Fetch all candidate activities before assay filtering, so excluded activity IDs
and reasons are auditable. Fetch distinct assay IDs once per invocation, in sorted
batches of at most 50 IDs; paginate every batch. Enrich only included molecules
and nonnull document IDs, also in batches of 50. Process eligible targets and
batch ID sets in lexicographic order, always using the fixed collection ordering
shown in the table. Skip empty ID batches. Single-resource detail endpoints
may be used for smoke checks; production should avoid one request per compound.
No independent UniProt call is needed for exact matching in this increment;
obsolete-accession resolution and reviewed-status verification remain out of scope.

### Pagination and joins

Require an object with the expected collection and `page_meta`. Require list items
to be objects. Pagination integers must be JSON integers, excluding booleans:
`offset >= 0`, `total_count >= 0`, `limit > 0`. Require `next` to be null or a
nonempty string; the initial returned offset must be 0. Subsequent returned offsets
must equal the requested next offset. The server may reduce the requested limit.

Resolve `next` against the current response URL. Permit only HTTPS, no credentials,
no fragment, host `www.ebi.ac.uk`, default/443 port, and the exact current collection
path. Require the original non-pagination parameters (including `order_by`) with
the same decoded values and no additional/duplicate parameters. Only limit and
offset may change; parse them as positive limit/nonnegative offset. The next offset
must equal current offset plus the number of rows returned and must increase.
Reject cycles, empty nonfinal pages, or more rows than the returned limit.

The first total_count is the expected number of distinct IDs for this collection
query (separately per target or batch). It must stay constant. A final `next: null`
page must bring the distinct-ID count to exactly total_count; earlier pages must
not exceed it. A structured empty collection with total_count 0 and next null is
valid. Identical repeated IDs may be collapsed, but they cannot conceal a missing
ID or justify a short final count. Conflicting full JSON objects for the same ID
fail validation (compare decoded objects, ignoring object-key order). Preserve
distinct activity IDs even when all other fields are identical. Database changes
or incomplete paging produce `failed`, error kind `schema`, with an explicit
incomplete-results message. There is no silent page/result cap.

Validate every target returned by accession lookup: a target with components that
do not contain the query is recorded as excluded, never included. Batch collections
must contain only requested IDs. Activities must reference the target requested by
that activity query, an assay ID, and a molecule ID. Resolve all distinct assay IDs
for candidate activities; resolve molecule/document IDs only for included activities.
A missing requested entity or a foreign/unexpected ID is a schema failure. An
optional document ID may be null. Use the activity document ID for measurement
references when present; otherwise fall back to the assay document ID. These IDs
may legitimately differ: ChEMBL permits measurements to have different sources
and documents from their associated assays (see the [ChEMBL data-model explanation](https://pmc.ncbi.nlm.nih.gov/articles/PMC6323927/)).
Do not treat this difference as a broken join or exclude the activity. Validate
both supplied IDs, but resolve/enrich the selected activity-first document only.
The original activity and assay document IDs remain traceable in their separate
raw responses, referenced by the request manifest and normalized activity/assay IDs.
Do not fetch activities for rejected targets. Audit counts therefore cover activities
of eligible targets, not activities of families/complexes rejected at resolution.

Use per-request three-attempt retry budgets and operation timeouts from the
BindingDB specification: 30 seconds per operation, retries for timeout/network
errors, 429, and 5xx; 1/2-second fallback backoff; Retry-After seconds or date;
stop rather than sleep above the 60-second limit. These are client policy choices,
not a claimed ChEMBL service quota. Throttle ChEMBL request starts to at most two
per second, sequentially, including retries and redirects. Space physical request
starts by at least 0.5 seconds using an injectable monotonic clock and sleeper; retry delay and throttle
are overlapping minimum delays, not additive sleeps.
Keep a run-local cache for repeated IDs; persistent caching is deferred.

Allow at most five HTTPS redirects within an attempt. Use the same authority
checks as pagination and restrict redirects to the same resource path and decoded
query (a trailing-slash normalization is permitted). Each retry restarts the logical
request URL and its redirect counter. Disable automatic HTTPX/transport retries.
Malformed JSON, HTML 200 responses, wrong types, missing envelopes, 404s, and
broken pagination fail the source. Empty bodies/JSON strings are schema failures;
only a valid completed collection can establish no matches. No scraping or dump
fallback. Preserve completed responses even when a later request fails.

Read `status.json` before and after the pipeline. If available versions disagree,
fail rather than mix releases. A missing version or failed status request emits
an explicit provenance warning and leaves `source_version` null; it does not
invalidate otherwise complete biological results. Status probes use the same
retry policy. A service reporting non-UP is a warning unless data requests fail.
Both status responses and warnings remain in the request manifest. Set
source_version only when both probes return the same nonempty version string.
Attempt the final probe after a required-request failure too, without replacing the original error;
on interruption propagate immediately. Two differing valid versions are a schema
failure even if the target collection was empty. Finish the source retrieved_at
timestamp after the final probe, using the existing injectable UTC clock format.
For a status response, require an object; missing/null/blank version is unavailable,
and a non-string version is invalid metadata. Either case produces warning code
`release_metadata_unavailable` with the probe role, just as a request/schema failure
on this optional endpoint does. A present non-UP status produces
`release_service_not_up`. Do not emit a no-match message for any failed source.

## Data model and traceability

### Validation and normalized measurements

Write saved-result schema version 2 for both sources. Keep the existing required
Measurement fields, permitting null `measurement_type` and `measurement_value` in
v2; all other existing types remain as specified for BindingDB. Append these fields
in exactly this order (also the appended TSV column order):

```text
source_record_id, target_id, target_accessions, assay_id, document_id,
assay_type, assay_description, target_organism, target_tax_id,
target_relationship, target_confidence, target_match_status,
variant_accession, variant_mutation, variant_isoform, activity_comment,
data_validity_comment, potential_duplicate, parent_compound_id,
measurement_relation, measurement_upper_value, measurement_text,
measurement_basis, source_values
```

All new fields are present in v2 JSON. Scalars default to null, target_accessions
to a new empty list, source_values to a new empty object, and target_match_status
to `unverified`. target_confidence is integer/null, potential_duplicate is
boolean/null (accept JSON bool or integer 0/1 only), and all other new scalars are
strings/null. BindingDB leaves new fields at defaults, including measurement_basis;
its existing value text, qualifiers, units, records and sorting remain unchanged.

ChEMBL validation rules:

- Entity IDs match `CHEMBL[0-9]+` with full ASCII matching. Activity IDs must be
  positive JSON integers, not booleans, and become decimal strings. Accessions
  must be nonempty strings when present; do not strip isoform suffixes to force
  identity. Required target_components is a list of objects; missing/wrong-type
  components or target_type fail schema validation. Empty components or null
  component accession/type make a target unverified and excluded. Component
  types and target types must otherwise be nonempty strings.
- Required activity fields are activity_id, assay_chembl_id, molecule_chembl_id,
  and target_chembl_id. Required assay fields are assay_chembl_id and
  target_chembl_id. Optional text may be missing/null/string; blank optional text
  becomes null, other strings are preserved. Wrong types fail, even on an activity
  later excluded. Unknown fields are ignored by normalization and retained raw.
- Missing/null assay confidence, relationship or type is valid unverified metadata
  and excludes its activities. Confidence must otherwise be an integer 0 through
  9, not bool. Unknown nonempty relationship/type codes are valid but excluded.
  If an activity supplies a nonnull assay_type, it must equal the assay's nonnull
  type; a missing activity type uses the assay. Null metadata is not a contradiction.
- `variant_sequence` may be absent/null/object. Map accession, mutation, isoform
  from that object. Fall back individually to activity assay_variant_accession and
  assay_variant_mutation when absent in the assay; two nonnull conflicting values
  fail schema validation. Normalize integer isoform to decimal string. Missing
  variant annotation is not a wild-type claim. Preserve the full variant object raw.
- Parse JSON decimals losslessly (Decimal, not binary float); reject NaN/Infinity.
  Numeric measurement fields accept finite JSON numbers or strings/null; normalize
  numbers with `str(Decimal)` (integer values with `str(int)`) and preserve strings.
  Numeric zero is present; never test measurement presence by truthiness. Objects,
  arrays, or booleans in numeric fields fail. Type, relation, units, text, descriptions,
  and names accept strings/null only. Document pubmed_id and tax IDs accept positive
  integers or ASCII digit strings/null, normalized to strings; isoform also permits 0.
- Null molecule_structures and molecule_hierarchy are valid; otherwise require
  objects with optional string structure fields/valid parent IDs. Missing document,
  name or structure metadata must not trigger an external-service fallback.

For included ChEMBL rows, source_record_id is activity_id; target_id is the resolved
target ID; target_accessions is `[query_accession]`; assay_id and document_id are
resolved IDs. Copy assay type/description/confidence/relationship, target organism
and tax_id, variant fields, activity comment/validity/duplicate flag and molecule
parent ID from their resolved entities. Set target_match_status to
`exact_accession_direct_assignment`; this describes source curation, not independent
experimental verification. A variant reference accession differing from the query
is annotation, not grounds to remap the target or include a different parent
accession; selection always uses the target component and assay assignment.

| Existing field | ChEMBL mapping |
| --- | --- |
| `compound_id` | `molecule_chembl_id`; do not collapse to parent/salt IDs |
| `compound_name` | Molecule `pref_name`, falling back to activity `molecule_pref_name`; otherwise null |
| `smiles`, `inchikey` | `molecule_structures.canonical_smiles`, `standard_inchi_key`; null structures allowed |
| `target_description` | Resolved target `pref_name` |
| `evidence_type` | `reported_binding_assay` for B; `reported_functional_assay` for F |
| `measurement_type/value/unit` | Coherent tuple chosen by the rules below |
| `pmid`, `doi` | Resolved document `pubmed_id` (decimal string), `doi`; optional |
| `source_url` | `https://www.ebi.ac.uk/chembl/web_components/explore/compound/CHEMBL_ID` |

Choose the standardized family if standard_type and standard_value are nonnull;
otherwise choose the published family if type and value are nonnull. A family
contains its own type, value, relation, units, upper_value, and text_value; units
or relation may be null and do not force fallback. Set measurement_basis to
`standard` or `published` accordingly and copy all selected family fields.
Never combine a published value with a standard unit. Retain relation and upper
bound separately; missing relation is unspecified, not implicit equality. For
text-only activity, measurement_type may use standard_type then type, value/unit
are null (as are relation and upper bound), measurement_basis is `qualitative`,
and text comes from standard_text_value then text_value then activity_comment.
If no type exists, display `Unspecified` while retaining null
in JSON. Preserve activity records even with no value/comment. Accordingly v2
permits null measurement_type/value. Render explicit missing values, never zero.
Do not convert IC50 to Ki, calculate pChEMBL, average, rank unlike types, or infer
missing units. Preserve source pchembl_value in source_values for audit only.
source_values is exactly an object with keys `published`, `standard`, `pchembl_value`: each of the
first two contains `type`, `value`, `relation`, `units`, `upper_value`, `text_value`
using normalized strings/null from that family, and pchembl_value is string/null.
Do not parse activity_comment as a compound ID or a number. Markdown and HTML
must display relations/ranges and qualitative text explicitly; preserve separate
machine-readable relation/value/unit fields in TSV.

### Result metadata and selection audit

Retain all existing SourceResult fields, adding `policy`, `warnings`, `counts`,
and `selection_audit`. policy is `exact-single-protein-v1` for ChEMBL and
`bindingdb-associated-v1` for BindingDB. warnings is a list of objects with string
`code` and `message`, and `role` (string/null, a manifest role for request warnings).
Errors retain the existing SourceError schema.

counts has exactly these keys: `complete` (bool), `candidate_targets`,
`eligible_targets`, `candidate_activities`, `included_activities`,
`excluded_activities`, `unverified_activities`, `unique_ligands` (integer/null).
On completed ChEMBL searches all counts are nonnegative integers; included plus
excluded equals candidate_activities. unverified_activities is the subset of
excluded activities with `unverified_target_assignment`, not an additional group.
On failure complete is false, candidate/eligible counts describe observations so
far, and included/excluded/unverified/unique counts are null; records is empty.
Do not present partially collected counts as complete.

selection_audit has exactly `outcome_reason` (string/null), `targets` and
`activities` (lists). Each target item has `target_id`, `included` (bool), and
`reasons` (sorted unique strings). Each activity item adds `activity_id`, `assay_id`
to those fields. Included entries have no reasons; excluded entries have all
applicable reasons. Sort targets by ID and activities by activity ID, lexicographically.
Raw responses retain all excluded values. Applicable target reasons are evaluated
independently; component_type is evaluated only for a sole nonnull component type.
Missing component metadata adds unverified_target_assignment instead of asserting
a known wrong component type. Stable reason codes are:

| Condition | Reason(s) |
| --- | --- |
| Non-SINGLE PROTEIN target | `target_type` |
| Component count is not one | `component_count` |
| Sole component type is non-PROTEIN | `component_type` |
| No component accession equals query | `component_accession_mismatch` |
| Component accession/type absent or component list empty | `unverified_target_assignment` |
| Assay relationship H | `homolog_assignment` |
| Other non-D nonnull relationship | `non_direct_assignment` |
| Nonnull confidence other than 9 | `confidence_not_9` |
| Nonnull assay type other than B/F | `unsupported_assay_type` |
| Null/missing assay relationship, confidence or assay type | `unverified_target_assignment` |

A mismatched target returned from accession lookup is excluded and audited. A
mismatched activity/assay join fails schema validation, not selection. Wrong data
types also fail rather than silently excluding. On success outcome_reason is null;
on no_matches choose `no_targets` if candidate_targets is zero, otherwise
`no_eligible_targets` if eligible_targets is zero, otherwise `no_activities` if
candidate_activities is zero, otherwise `all_activities_excluded`. On failure it is `incomplete`; any audit
entries already evaluated remain available, and absence from a partial audit is
not a selection decision.

Target audit size equals candidate_targets; included target entries equal
eligible_targets. Completed activity audit size equals candidate_activities;
true/false included flags give included/excluded counts. Derive unique_ligands
from included records; each record's source_record_id occurs once. Per-reason
report counts count entries containing that reason, so their sum may exceed total
excluded activities. Do not sum reason counts as a total.

For BindingDB, candidate_targets/eligible_targets/unverified_activities are null,
candidate_activities equals len(records), included_activities equals len(records),
excluded_activities is 0, and unique_ligands is its distinct monomer count when
completed. On failure use the failure/null rules above. Its audit lists are empty;
outcome_reason is null on success, `no_activities` on no_matches, `incomplete` on
failure. Do not fabricate per-row decisions or IDs for BindingDB.

### Raw-response and request manifest contract

Keep BindingDB's existing final raw_body and raw filename behavior. ChEMBL returns
raw_body null and an ordered list of HTTP exchanges, retained in memory until
storage completes. An exchange is one physical HTTP request; each redirect hop
is a separate exchange;
a transport failure is an exchange with no response. Logical requests get increasing
request_id values starting at 1; retries share that ID, with attempt 1 through 3 and redirect_hop
0 through 5. Every exchange gets a monotonically increasing sequence starting at 1.

`sources/chembl-requests.json` contains exactly `schema_version: 1`,
`accession` (normalized query), `policy`, `warnings`, and `exchanges`. Each exchange has exactly:

```text
sequence (int), request_id (int), role (string), requested_ids (list[string]),
attempt (int), redirect_hop (int), url (full encoded request URL),
parameters (object[string,string]), retrieved_at (UTC string),
http_status (int|null), transport_error (string|null),
location (string|null), raw_path (string|null), sha256 (string|null)
```

Roles are `status_before`, `targets`, `activities`, `assays`, `molecules`,
`documents`, `status_after`. requested_ids is sorted; it contains the query accession
for targets, target ID for activities, batch entity IDs for joins, and is empty
for status. parameters is the actual decoded query, including offset/limit when
present. Store only a concise transport error, never a traceback/body in error text.

Each complete body is stored verbatim as bytes at
`chembl-raw/NNNNNN.json` if strict JSON parses, otherwise `chembl-raw/NNNNNN.txt`,
where NNNNNN is the exchange sequence padded to at least six digits. raw_path is
relative to `sources/`; sha256 hashes those exact bytes. Retain complete empty,
error, and redirect bodies too. No response means null raw_path/sha256/http_status.
This distinguishes a 200 empty body from transport failure. Do not serialize bytes
into SourceResult or the manifest; storage pairs exchange metadata with body bytes.

ChEMBL SourceResult.request is exactly `{"url": BASE + "target.json",
"parameters": {"target_components__accession": ACCESSION, "limit": "100",
"offset": "0", "order_by": "target_chembl_id"},
"manifest": "chembl-requests.json"}`; manifest is relative to sources/.
BindingDB's existing request object is unchanged. SourceResult.source is `ChEMBL`
or `BindingDB`. JSON is UTF-8, indent 2, newline terminated; write decimal-derived
strings, never JSON Decimal objects. source_version follows the status rules
above. Each warning is copied to the source result and manifest; order warnings by probe order and then warning code,
deduplicating identical role/code/message triples.

## Source states, orchestration, and outputs

ChEMBL returns `success` when at least one activity survives selection;
`no_matches` after a complete successful pipeline with none, including no exact
targets, no activities, or all candidates excluded. Record which case occurred
and report excluded counts. A required request/validation failure produces
`failed`, records empty, and explicit `transport`, `http`, or `schema` error.
Retain candidate evidence/manifest for diagnosis but do not publish a partial
ligand count. Enrichment failures are required-request failures in this first
increment; missing optional fields are not. Catch expected service failures as
data, not unexpected programming errors.

Change orchestration to a run result containing ordered source outcomes
(BindingDB, ChEMBL). Keep each source independent; do not abandon the second when
the first fails, and never reuse cached results from an earlier invocation.
Exit 0 if all sources completed (including no_matches), 1 if any source or output
publication failed, 2 for CLI errors/explore, 130 for interruption. On a source
failure still save and display the other source's current results.

Output layout:

```text
results/ACCESSION/
  report.md
  ligands.tsv
  bindingdb-ligands.html
  chembl-ligands.html
  sources/
    bindingdb.json
    bindingdb.raw.json       # existing raw variant rules retained
    chembl.json
    chembl-requests.json
    chembl-raw/              # numbered response bodies referenced by manifest
```

`report.md` begins with a consolidated summary and this exact table:

```markdown
| Database | Ligands Count | Ligands Report |
| --- | --- | --- |
| BindingDB | N | [View ligands](bindingdb-ligands.html) |
| ChEMBL | M | [View ligands](chembl-ligands.html) |
```

Count distinct compound IDs independently per source after that source's policy.
A failed source shows `Unavailable`, not 0; a successful no-match source shows 0.
Show source status, release, query accession, selection policy, warnings, and
exclusions alongside the table. State that ChEMBL uses stricter assay-target
selection than BindingDB, so counts are not directly equivalent. Do not add
source counts into a claimed globally unique molecule count: identifiers and
possibly underlying measurements overlap across databases.

Generate one local HTML file per source with one row per source compound ID,
a visible SMILES column, expandable complete measurements, supplied names, target/assay links, relation,
value/range or qualitative text, units, organism, target assignment, variant,
quality flags, and references. Compound URLs use the mapping above; target, assay and document
URLs use `https://www.ebi.ac.uk/chembl/web_components/explore/` plus
`target/ID`, `assay/ID`, or `document/ID`. No such links are fabricated when
BindingDB lacks those identifiers. Keep the existing digits-only PMID and
`10.` DOI checks and identifier URL encoding. Escape all text and attribute values;
construct links from validated IDs on official hosts. Generalize the existing HTML renderer: it currently
hardcodes BindingDB compound links. No external assets or scripts are required.
BindingDB retains its unverified-evidence notice; ChEMBL states "Exact accession
and direct assay-target assignment; physical binding not independently verified."
Always state that the ChEMBL count includes inactive, inconclusive, qualitative,
and functional assay records. Display warnings, selection policy, no-match reason,
and failed-source counts as Unavailable in HTML too. Show rejected/unverified
activity counts and link to `sources/chembl.json` for the full audit; do not put
excluded compounds in the counted ligand table.
The SMILES column contains each distinct nonblank supplied string for that
compound, sorted lexicographically, once per ligand row. Preserve strings exactly,
HTML-escape them and render each in a separate selectable monospace code block.
Wrap long strings visually without inserting characters or clipping the string;
use `user-select: all` for convenient copying. Missing structures display
`Not supplied`. Do not join multiple SMILES into one chemical string.
Both source HTML files link back to `report.md` and the combined TSV.

End every source section in terminal output with one empty line, separating
successive databases and the final Results/output-path summary.
Terminal output shows the accession, followed by each source's status, unique
ligand count or Unavailable, measurement count or Unavailable, and at most five
reported targets and five measurement types by record count. Null type is
`Unspecified` before grouping/sorting. Retain existing label sanitization, ordering,
truncation and omitted-label counts. Show no-match reasons, warnings, failures,
and source-specific evidence notices without individual measurements or full audit
entries. End with the output folder and both HTML filenames. Log source failures
to stderr; print summary only after publication succeeds.

Retain the full measurement tables in Markdown for compatibility, separated by
source, but put the compact database table and HTML links first. Keep existing
TSV columns in their current order and append v2 fields in the order listed above;
encode target_accessions and source_values as compact JSON cells (ensure_ascii=False,
sort_keys=True for objects); encode scalar null as empty, bool as `true`/`false`,
and confidence as decimal digits. Keep UTF-8, csv quoting and line endings from
the existing specification. The TSV remains header-only if both sources have no
records; no synthetic failure rows. Repeat source, query accession, and source retrieval time per row. Source order is BindingDB then
ChEMBL. Sort ChEMBL rows by compound ID, target ID, assay ID, and activity ID
lexicographically. Preserve BindingDB's existing stable sort.

Generalize staging/publication for both sources. Render and stage all artifacts
before replacing generated destinations. Publish in this order: BindingDB raw body
and ChEMBL raw exchanges; ChEMBL manifest; BindingDB then ChEMBL HTML; combined TSV;
combined Markdown; BindingDB then ChEMBL source JSON. After publishing current files,
remove obsolete BindingDB raw variants and ChEMBL raw files that match
`[0-9]{6,}\.(json|txt)` but are not referenced by the new manifest. Never recursively
remove the raw directory or touch other filenames. An empty or failed rerun must replace old ChEMBL HTML/JSON and remove stale counted rows from
combined outputs. Preserve other files. Atomicity remains per file, not per run;
report output failures as potentially mixed results and exit 1. Concurrent writes
to the same accession directory remain unsupported.

## Deterministic acceptance example

Use synthetic fixtures with query P00533, one eligible single-protein target, and
one excluded two-component complex containing P00533. Return these five distinct
activities for the eligible target; all joins resolve and no metadata contradicts:

| Activity | Molecule | Assay metadata | Expected selection |
| --- | --- | --- | --- |
| 1 | CHEMBL1 | B, D, 9, quantitative | Include |
| 2 | CHEMBL1 | F, D, 9, L858R annotation | Include, same unique molecule as activity 1 |
| 3 | CHEMBL2 | B, H, 8 | Exclude: homolog_assignment, confidence_not_9 |
| 4 | CHEMBL3 | B, null relationship, 9 | Exclude: unverified_target_assignment |
| 5 | CHEMBL4 | B, D, 9, comment Inactive and null numeric values | Include qualitative record; do not claim binding |

Expected counts: complete true, candidate_targets 2, eligible_targets 1,
candidate_activities 5, included_activities 3, excluded_activities 2,
unverified_activities 1, unique_ligands 2; status success; outcome_reason null.
The ChEMBL database table says 2, HTML has two ligand rows and three expanded
activity entries, and combined TSV has three ChEMBL rows. Complex activities
are not fetched. The target rejection reasons are component_count and target_type.
This fixture is synthetic; its molecule/activity numbers are not live data claims.

## Implementation sequence and required tests

Implementation interfaces (dataclasses and ordinary functions; no model-validation
or ChEMBL SDK dependency):

```text
fetch_bindingdb(accession, *, client, sleep, now) -> SearchOutcome
fetch_chembl(accession, *, client, sleep, now, monotonic) -> SearchOutcome
run_simple(accession, *, client, sleep, now, monotonic=time.monotonic) -> RunOutcome
render_summary(run: RunOutcome, *, folder: Path) -> str
render_report(run: RunOutcome) -> str
render_tsv(run: RunOutcome) -> str
render_html(result: SourceResult) -> str
save_results(run: RunOutcome, *, root: Path) -> Path
```

RunOutcome has accession and ordered sources (list[SearchOutcome]). SearchOutcome
keeps result and raw_body, adding exchanges (list[HttpExchange], default empty).
HttpExchange holds the manifest fields plus nullable response body bytes; storage
computes raw_path and sha256. fetch_bindingdb returns its adapter-level schema-1
result; run_simple upgrades it to v2 metadata/defaults without modifying the adapter.
No other adapter needs ChEMBL HTTP/provenance logic. CLI supplies time.monotonic,
time.sleep and the UTC clock; adapters do not close caller-owned clients.

Add `source_result_from_dict(data) -> SourceResult` for v1/v2 saved JSON. For v1,
require existing keys/types, fill v2 measurement defaults and derive BindingDB
metadata as above; for v2 require the defined fields/types. Reject unknown schema
versions. The returned in-memory result uses schema_version 2. The v1 reader
accepts only the existing BindingDB source contract; it must not guess a ChEMBL
policy for a version-1 file. Preserve original version-1 request/provenance fields.
This is a library/fixture reader, not a new CLI input or permission to
mix old runs. Do not rewrite existing files unless explicitly saving a new run. The source renderers
and storage also accept single SourceResult/SearchOutcome values where applicable
for library compatibility; the CLI always supplies a RunOutcome containing both sources.

Keep network access and parsing separate enough for individual fixture tests;
renderers do not fetch data, and the search service does not publish files.
Source-wide expected failures return result data; only interruption/programming
errors propagate. Update relevant existing BindingDB specification interfaces and
README in the implementation change to describe this contract.

The acceptance matrix below is the regression contract for the implemented integration.

1. Introduce v2 models, collection validation, source-aware URLs, and multi-response
   provenance without changing BindingDB parsing/count policy.
2. Add `chembl.py` with accession resolution, paginated retrieval, assay selection,
   batched molecule/document enrichment, retry/throttle handling, and raw manifest.
3. Generalize `search.py`, `storage.py`, and report renderers for multiple outcomes.
   Update CLI status aggregation and README coverage/output examples.
4. Add mocked HTTPX tests and small fixtures; run relevant tests and the full
   pytest suite before completion. Live checks are optional and never assert
   immutable live counts. Documentation-only revisions do not require pytest.

Acceptance cases must cover:

- Exact request filters (including verified `__in` and `order_by` parameters),
  all targets/pages consumed, and multiple eligible IDs; no targets versus targets with zero activities versus all activities excluded.
- P00533-style lookup containing family/chimera/interaction results: only exact
  single-protein targets eligible. P0DTC1 no coverage must not query P0DTD1.
- Exact target ID but assay H/8 excluded (activity 32260 pattern); D/9 B and F
  included; null/unknown assignments excluded and counted; malformed metadata
  fails rather than silently excluding; assay/activity target mismatch fails.
- Different activity/assay documents with activity-first references; null activity
  document falls back to the assay; null documents remain valid; invalid document
  IDs and missing selected documents fail.
- Variants and null variant metadata; qualitative/inactive/unknown endpoints;
  source qualifiers, ranges, zero values, null values, coherent units, and
  source flags preserved. Two distinct activities for one ligand count once.
- Multiple batches and pages for all collections; batch boundary at 50 IDs;
  missing referenced records; repeated IDs with identical/conflicting data;
  malicious next URLs, cycles, changing totals, and failures after earlier pages.
- Null structures/names/documents; real nullable fields versus wrong types;
  integer PMID conversion; no accidental parent/salt deduplication.
- Timeout/429/5xx retries, Retry-After, redirects, maximum delay, throttling,
  malformed bodies, empty strings, per-request retry limits, and version-change
  failure. Status metadata outage yields warning/unknown version, not no matches.
- One source fails while the other succeeds; all fail; all complete with no matches;
  exit codes; no stale results; current-run provenance and every raw page saved.
- Exact Markdown table/relative links, safe source-specific compound URLs,
  readable HTML, escaped untrusted values, source-separated counts, and TSV v2
  round trips. BindingDB rendering/units/counts and v1/v2 reader defaults and
  unknown-version regressions.
- Publication/staging failures, obsolete raw cleanup, success-to-empty/failed
  reruns, interruption cleanup, and preservation of unrelated files.

Completion means both online sources work through the installed CLI, each has a
traceable source result and readable local report, filters are auditable, all
required tests pass, and the implementation specifications/README match the code.


## Implementation verification

A live installed-command smoke check on November 4, 2026 used `ligandscout P00325`
and completed with exit code 0. BindingDB returned 24 distinct compounds and 24
measurements. ChEMBL_37 returned 24 candidate activities; the direct-assignment
policy included one activity/compound and excluded 23 activity records
with both homolog_assignment and confidence_not_9. These are observations, not
fixed regression expectations. Both source HTML files and consolidated table
links resolved locally, both saved source JSON files loaded through the v2 reader,
and every ChEMBL raw response matched its manifest SHA-256.
