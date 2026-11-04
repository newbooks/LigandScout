# LigandScout

Find known and potential ligands for a protein using its **Swiss-Prot accession ID** and available ligand resources.

- **Simple mode:** Collect reported ligands and supporting evidence.
- **Explore mode:** Infer potential ligands from PClass embedding neighbors and PClass-Dark cluster members.

## Ligand Databases

Planned resources include ligand, annotation, and compound databases. API and download options were checked in November 2026.

| Database | Coverage and use | Online API | Downloads |
| --- | --- | --- | --- |
| BindingDB | Protein–small-molecule affinity measurements | [REST](https://www.bindingdb.org/rwd/bind/BindingDBRESTfulAPI.jsp) | [TSV, SDF, database dumps](https://www.bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp) |
| ChEMBL | Compound–target bioactivity and assays | [REST](https://chembl.gitbook.io/chembl-interface-documentation/web-services/chembl-data-web-services) | [Database dumps, SDF](https://chembl.gitbook.io/chembl-interface-documentation/downloads) |
| IUPHAR/BPS Guide to Pharmacology | Curated ligands and pharmacological interactions | [REST; API key required](https://www.guidetopharmacology.org/webServices.jsp) | [Ligand and interaction tables](https://www.guidetopharmacology.org/download.jsp) |
| PubChem / BioAssay | Chemical records and screening results | [PUG REST](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest) | [Bulk data](https://ftp.ncbi.nlm.nih.gov/pubchem/) |
| RCSB PDB | Ligands observed in experimental structures | [Search / Data APIs](https://www.rcsb.org/docs/programmatic-access/web-apis-overview) | [Structures and ligand files](https://www.rcsb.org/docs/programmatic-access/file-download-services) |
| BRENDA | Enzyme substrates, products, inhibitors, activators, and kinetics | [SOAP; registration required](https://www.brenda.dsmz.de/soap.php) | [Text / JSON; accept license](https://www.brenda-enzymes.de/download.php) |
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

| Database | Ligands |
| --- | --- |
| ... | ... |

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
