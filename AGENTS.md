# Repository Ground Rules

## Branches and approval

- Never edit files or create commits on `main`. Switch to a working branch before making changes.
- On branches other than `main`, proceed with task-related code changes, documentation updates, and tests without requesting additional approval. This standing authorization remains subject to tool permissions. It does not authorize destructive actions, merges into `main`, publishing, or external communications.
- Preserve unrelated user changes.

## Ligand search: API access

- Limit source integrations to BindingDB, UniProtKB/Swiss-Prot, and RCSB PDB. BindingDB and UniProt are implemented. RCSB CCD structure enrichment is implemented; experimental complex/contact search remains planned.

- Always implement ligand searches through the resources' online APIs. This is the agreed architecture for occasional batches of about 20 accessions.
- Do not implement searches against downloaded database dumps or local database mirrors unless the user explicitly changes this decision. Local caching of API responses is allowed.
- Start from a Swiss-Prot accession and resolve any additional identifiers required by each API through API calls.
- Handle pagination, rate limits, timeouts, and retries. Distinguish a successful search with no matches from a failed or unavailable service.
- Retain source identifiers, query parameters, retrieval dates, and database versions when available so results remain traceable.

## Implementation specifications

- Follow the [CLI and BindingDB implementation specification](docs/specs/cli-bindingdb.md) for the implemented simple-mode command. The UniProt extension is defined in the [UniProt integration specification](docs/specs/uniprot-integration.md).
- Whenever an algorithm changes, update the relevant implementation specification as part of the same change without waiting for a separate request. Keep the documented inputs, outputs, steps, assumptions, and limitations consistent with the code.
- If no implementation specification exists, create one under `docs/specs/` for the affected functionality.

## Tests

- Every code change must include new or updated pytest tests for the affected functionality. Cover meaningful edge cases and add regression tests for bug fixes.
- Run the relevant pytest tests whenever you add or modify code. Resolve failures caused by the change before considering the work complete.
- Report which tests ran and their results. If tests cannot run, explain the blocker and do not claim they passed.
- Documentation-only changes do not require pytest tests.
