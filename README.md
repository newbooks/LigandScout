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

Simple mode queries **BindingDB and ChEMBL** through their online APIs. It accepts one canonical UniProt accession, assumed to belong to Swiss-Prot. Running without an accession displays help. `--explore` is reserved and is not implemented yet.

Results are saved relative to the current directory:

```text
results/P00699/
├── sources/
│   ├── bindingdb.json
│   ├── bindingdb.raw.json
│   ├── chembl.json
│   ├── chembl-requests.json
│   └── chembl-raw/
├── bindingdb-ligands.html
├── chembl-ligands.html
├── ligands.tsv
└── report.md
```

Open `report.md` for a database/count/report table, then follow each source's local HTML link. HTML lists one row per unique source compound ID, with a visible, selectable SMILES column, expandable measurements and links to database entries. `ligands.tsv` contains all included measurements from both sources. Counts are separate per source: the same compound or experiment can appear in both databases, so the counts are not a globally unique total.

ChEMBL requires an exact UniProt component accession on a single-protein target and a direct assay-target assignment with confidence 9. It excludes homolog assignments, families and complexes. Both binding and functional assays are included, including inactive, inconclusive and qualitative records; these are tested compounds, not necessarily proven binders. Variants are retained and labeled where annotated. A query with no exact ChEMBL coverage does not silently search a related accession. BindingDB retains its API associations without independently verified target identities. IC50 alone does not establish physical binding.

Source JSON uses schema version 2 and records status, retrieval time, policy and provenance. ChEMBL additionally saves selection decisions, excluded/unverified activity counts, release information when available, and a manifest linking every HTTP exchange to its raw response and SHA-256. BindingDB preserves its final raw response (`.txt` for non-JSON bodies). The TSV keeps the original columns and appends assay, target, variant, selection and original-value metadata. Missing units remain unknown.

Rerunning an accession makes fresh API requests and **overwrites generated results**, replacing stale rows even after an empty or failed search. Sources run independently: one failure still allows the other source's results to be saved. Failed sources show `Unavailable`, distinct from a successful zero-match search. Ambiguous BindingDB empty bodies are retried; ChEMBL requires a complete structured collection. Release metadata outages produce warnings without invalidating otherwise complete ChEMBL results.

Use `--log-level DEBUG` for detailed logs. Exit codes are 0 when both sources complete (including no matches), 1 if any source or output fails, 2 for invalid input or the explore placeholder, and 130 for interruption. Publication is atomic per file, not across the output directory; if saving fails or is interrupted, rerun because the folder may contain mixed results. Concurrent searches into the same accession folder are unsupported. Annotation, other sources, and explore-mode features below remain planned.

Implementation contracts: [BindingDB](docs/specs/cli-bindingdb.md) and [ChEMBL / multi-source reporting](docs/specs/chembl-integration.md).

## Ligand Databases

BindingDB and ChEMBL are implemented. Additional planned resources include ligand, annotation, and compound databases. API and download options were checked in November 2026.

| Database | Coverage and use | Online API | Downloads |
| --- | --- | --- | --- |
| BindingDB | Protein–small-molecule affinity measurements | [REST](https://www.bindingdb.org/rwd/bind/BindingDBRESTfulAPI.jsp) | [TSV, SDF, database dumps](https://www.bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp) |
| ChEMBL | Compound–target bioactivity and assays | [REST](https://chembl.gitbook.io/chembl-interface-documentation/web-services/chembl-data-web-services) | [Database dumps, SDF](https://chembl.gitbook.io/chembl-interface-documentation/downloads) |
| PubChem / BioAssay | Chemical records and screening results | [PUG REST](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest) | [Bulk data](https://ftp.ncbi.nlm.nih.gov/pubchem/) |
| RCSB PDB | Ligands observed in experimental structures | [Search / Data APIs](https://www.rcsb.org/docs/programmatic-access/web-apis-overview) | [Structures and ligand files](https://www.rcsb.org/docs/programmatic-access/file-download-services) |
| UniProtKB/Swiss-Prot | Binding sites, cofactors, target cross-references, and subunit structure | [REST](https://www.uniprot.org/help/programmatic_access) | [Release files](https://www.uniprot.org/help/downloads) |
| UniProt-GOA | GO molecular functions and relations describing how a protein contributes to an activity | [QuickGO](https://www.ebi.ac.uk/QuickGO/api/index.html) | [GOA annotation files, including GAF](https://www.ebi.ac.uk/GOA/downloads) |
| Rhea | Reaction substrates and products, linked through protein annotations | [REST / SPARQL](https://www.rhea-db.org/help/api) | [Reaction data](https://www.rhea-db.org/help/download) |
| ChEBI | Compound identifiers, structures, and roles for normalization | [REST](https://www.ebi.ac.uk/chebi/tools) | [SDF, ontology, database dumps](https://www.ebi.ac.uk/chebi/downloads) |

Assay hits, reaction participants, and structural observations require interpretation before they can support a claim of direct binding. Record unavailable sources and retrieval dates. Inferred candidates require validation.

### Functional context: alone or in a complex

Use Swiss-Prot **Subunit structure** (`CC   -!- SUBUNIT:` in the flat file), **Function**, and their supporting evidence to assess whether the query protein performs an activity as a monomer, a homomer (multiple copies of the same protein), or a heteromeric complex (different proteins). Report this context for each activity where evidence is available; otherwise report it as unknown.

Supplement this assessment with UniProt-GOA molecular-function annotations. In [GAF 2.2](https://geneontology.org/docs/go-annotation-file-gaf-format-2.2/), column 4 explicitly records `enables` or `contributes_to`, optionally prefixed with `NOT|`. `contributes_to` can capture a subunit's contribution to an activity of a complex; **`enables` does not establish that the protein functions alone**. Preserve negation and annotation evidence rather than treating either relation as a definitive oligomeric-state label.

Swiss-Prot [flat-file GO cross-references](https://web.expasy.org/docs/userman.html) (`DR   GO; ...`) contain the GO identifier, aspect/term, and evidence/source, but do not preserve these relations as structured fields. Obtain the relations from GOA and match annotations by UniProt accession, retaining isoform information when present. Check the `!gaf-version:` header: in [GAF 2.1](https://geneontology.org/page/go-annotation-file-gaf-format-21), column 4 is an optional qualifier field that supports `contributes_to`, while the ordinary positive molecular-function relation is implicit rather than an explicit `enables` value.

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
