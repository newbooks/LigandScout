# LigandScout

Find known and potential ligands for a protein using its **Swiss-Prot accession ID** and available ligand resources.

- **Simple mode:** Collect reported ligands and supporting evidence.
- **Explore mode:** Infer potential ligands from PClass embedding neighbors and PClass-Dark cluster members.

## Installation

```bash
pip install -e .
ligandscout --help
```

To install test dependencies, use `pip install -e '.[test]'`.

## Search

```bash
ligandscout P00699
ligandscout P00176 --simple
```

The selected resources are **BindingDB, UniProtKB/Swiss-Prot, and RCSB PDB**. Simple mode queries **BindingDB and UniProt** independently. UniProt ligands are enriched with SMILES from RCSB chemical components only when chemical identifiers match exactly; a separate experimental PDB complex search remains planned. It accepts one canonical UniProt accession, assumed to belong to Swiss-Prot. Running without an accession displays help. `--explore` remains a placeholder.

Results are saved relative to the current directory:

```text
results/P00699/
├── sources/
│   ├── bindingdb.json
│   ├── bindingdb.raw.json
│   ├── uniprot.json
│   ├── uniprot-requests.json
│   └── uniprot-raw/
├── index.html
├── bindingdb-ligands.html
├── uniprot-ligands.html
├── ligands.tsv
└── report.md
```

The terminal displays a compact summary. Open `results/<accession>/index.html` for the summary, source status, ligand and measurement counts, retrieval times, and links to each database’s results. Individual database pages link back to this index. `report.md` remains available as a detailed Markdown report. BindingDB HTML lists unique compounds with selectable SMILES and expandable measurements. UniProt HTML shows annotated ligands/cofactors, binding-site coordinates and residue sequences, evidence, and available chemical structures. Generic sites and reaction/function context appear separately when available. The UniProt page presents binding sites as bullet lists in the ligand table and short source-text extracts with links to the original UniProt entry; full JSON stays in the saved data files. `ligands.tsv` combines the sources, keeping affinity measurements separate from annotation records; site and structure details are JSON-valued columns.

Results are API-associated records, not confirmed binding partners. BindingDB target identities and assay evidence have not been independently verified; IC50 alone does not establish physical binding. Missing units remain unknown. UniProt annotations may be experimental or inferred; their evidence remains attached to the original site or cofactor. Binding sites use 1-based UniProt sequence coordinates. RCSB CCD supplies chemical structures only, not evidence of binding to this target. SMILES can be unavailable or ambiguous, especially for different charge states; the report explains the reason. Reaction participants are not automatically labeled as binders.

Source JSON uses schema version 4; the reader also supports BindingDB v1/v2/v3 files. Saved results retain retrieval time, query provenance, policy, status, and available measurement/target metadata. BindingDB preserves its final raw response (`.txt` for non-JSON bodies). UniProt preserves all entry/CCD response bodies with request provenance and SHA-256 hashes in its request manifest.

Rerunning an accession makes fresh API requests and **overwrites generated source results and combined reports**, replacing stale rows even after an empty or failed search. Failed searches show `Unavailable`, distinct from a successful zero-match search. Ambiguous empty bodies are retried and recorded as failures if they persist. Previously saved files from removed integrations are left untouched but are never read into new results or linked from new reports.

A source failure still allows the other source to complete and be saved. CCD enrichment outages retain UniProt annotations with warnings and missing structures. Use `--log-level DEBUG` for detailed logs. Exit codes are 0 for completed searches (including no matches), 1 for service/output failures, 2 for invalid input or the explore placeholder, and 130 for interruption. Publication is atomic per file, not across the directory; rerun after an output failure or interruption. Concurrent searches into the same accession folder are unsupported.

Implementation contracts: [CLI and BindingDB](docs/specs/cli-bindingdb.md), [UniProt ligands and sites](docs/specs/uniprot-integration.md).

## Ligand Databases

Only these three resources are in scope. All integrations use online APIs.

| Database | Purpose | Status | API |
| --- | --- | --- | --- |
| BindingDB | Protein–small-molecule activity and affinity measurements | Implemented | [REST](https://www.bindingdb.org/rwd/bind/BindingDBRESTfulAPI.jsp) |
| UniProtKB/Swiss-Prot | Protein identity, processed chains, biological ligands, cofactors, binding sites, and subunit context | Implemented | [REST](https://www.uniprot.org/help/programmatic_access) |
| RCSB PDB | Experimental protein–ligand complexes and binding-site context | CCD structure enrichment implemented; complex search planned | [Search / Data APIs](https://www.rcsb.org/docs/programmatic-access/web-apis-overview) |

Assay hits, reaction participants, and structural observations require interpretation before they can support a claim of direct binding. Record unavailable sources and retrieval dates. Inferred candidates require validation.

### Functional context: alone or in a complex

Use Swiss-Prot **Subunit structure** (`CC   -!- SUBUNIT:` in the flat file), **Function**, and their supporting evidence to assess whether the query protein performs an activity as a monomer, a homomer (multiple copies of the same protein), or a heteromeric complex (different proteins). Report this context for each activity where evidence is available; otherwise report it as unknown.

UniProt annotations alone may not resolve every activity-specific complex relationship. Report uncertainty explicitly rather than inferring that a protein acts alone from a GO molecular-function term.

## Planned outputs

Produce a human-readable report with JSON and TSV exports. Include source links and supporting evidence, and distinguish reported ligands from inferred candidates.

### Simple mode report example

- Query accession: ...
- Description: ...
- GO molecular function terms: ...
- Functional context: monomer / homomer / heteromeric complex / unknown, with activity-specific GO relations, partner or complex names, and supporting evidence where available
- Sequence length: ...

| Database | Ligands Count | Ligands Report |
| --- | --- | --- |
| ... | ... | ... |

### Explore mode report example

Add two candidate tables to the simple-mode report. Include the number of neighbors or members with usable data alongside each percentage.

Potential ligands inferred from the top 10 embedding neighbors (query side):

| Database | Ligand | Reported percentage |
| --- | --- | --- |
| ... | ... | ... |

Potential ligands inferred from cluster members (target side):

| Database | Ligand | Reported percentage |
| --- | --- | --- |
| ... | ... | ... |

## License

See [LICENSE](LICENSE) for the repository license. External data and tools remain subject to their own terms.
