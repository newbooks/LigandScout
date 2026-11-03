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
| UniProtKB/Swiss-Prot | Binding sites, cofactors, and target cross-references | [REST](https://www.uniprot.org/help/programmatic_access) | [Release files](https://www.uniprot.org/help/downloads) |
| Rhea | Reaction substrates and products, linked through protein annotations | [REST / SPARQL](https://www.rhea-db.org/help/api) | [Reaction data](https://www.rhea-db.org/help/download) |
| ChEBI | Compound identifiers, structures, and roles for normalization | [REST](https://www.ebi.ac.uk/chebi/tools) | [SDF, ontology, database dumps](https://www.ebi.ac.uk/chebi/downloads) |

Assay hits, reaction participants, and structural observations require interpretation before they can support a claim of direct binding. Record unavailable sources and retrieval dates. Inferred candidates require validation.

## Planned outputs

Produce a human-readable report with JSON and TSV exports. Include source links and supporting evidence, and distinguish reported ligands from inferred candidates.

### Simple mode report example

- Query accession: ...
- Description: ...
- GO molecular function terms: ...
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
