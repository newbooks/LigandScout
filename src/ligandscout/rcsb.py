"""Exact UniProt entity selection and RCSB-reported ligand contacts."""

import json
import logging
import math
import re
from urllib.parse import quote

from ligandscout.accessions import normalize_accession
from ligandscout.api_http import Audit, request_json, timestamp
from ligandscout.models import Measurement, SearchOutcome, SourceError, SourceResult
from ligandscout.rcsb_schema import counts, validate_result

SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
CORE = "https://data.rcsb.org/rest/v1/core/"
LOG = logging.getLogger(__name__)
LIMITS = [
    "Experimental PDB nonpolymer entities only; polymer partners, branched carbohydrates and polymer modifications are excluded.",
    "Contacts are RCSB-reported neighborhoods; assemblies, symmetry mates, conformers and models may not be fully covered.",
    "Structural proximity does not establish affinity or physiological binding.",
]
AMINO = dict(
    zip(
        [
            "ALA",
            "ARG",
            "ASN",
            "ASP",
            "CYS",
            "GLN",
            "GLU",
            "GLY",
            "HIS",
            "ILE",
            "LEU",
            "LYS",
            "MET",
            "PHE",
            "PRO",
            "SER",
            "THR",
            "TRP",
            "TYR",
            "VAL",
        ],
        "ARNDCQEGHILKMFPSTWYV",
        strict=True,
    )
)


class RCSBError(ValueError):
    """Expected schema/HTTP failure with source provenance."""

    def __init__(self, message, error=None):
        super().__init__(message)
        self.error = error


def demand(condition, message):
    if not condition:
        raise RCSBError(message)


def string(value):
    return isinstance(value, str) and bool(value.strip())


def positive(value):
    return type(value) is int and value > 0


def array(obj, key, *, empty=True):
    value = obj.get(key, [])
    demand(isinstance(value, list) and (empty or bool(value)), f"Invalid {key} list.")
    return value


def identifiers(obj, key, *, empty=True):
    values = array(obj, key, empty=empty)
    demand(
        all(string(x) for x in values) and len(set(values)) == len(values),
        f"Invalid {key} identifiers.",
    )
    return values


def mapping(entity, position, accepted, preferred, sequence_entry):
    root = entity.get("rcsb_polymer_entity_align", [])
    base = {
        "uniprot_accession": None,
        "uniprot_position": None,
        "mapping_reason": "no_alignment",
        "uniprot_residue": None,
        "sequence_check": "not_checked",
        "alignment_segments": [],
    }
    try:
        demand(isinstance(root, list), "Invalid alignment collection.")
        matched = []
        for item in root:
            if not isinstance(item, dict):
                raise RCSBError("Invalid alignment object.")
            if (
                item.get("reference_database_name") != "UniProt"
                or item.get("reference_database_accession") not in accepted
            ):
                continue
            segments = array(item, "aligned_regions", empty=False)
            for seg in segments:
                demand(
                    isinstance(seg, dict)
                    and all(
                        positive(seg.get(k))
                        for k in ("entity_beg_seq_id", "ref_beg_seq_id", "length")
                    ),
                    "Invalid alignment segment.",
                )
                demand(
                    seg["entity_beg_seq_id"] + seg["length"] - 1
                    <= entity["entity_poly"]["rcsb_sample_sequence_length"],
                    "Alignment outside entity.",
                )
            matched.append(item)
        if not matched:
            return base
        for accession in dict.fromkeys((preferred, *accepted)):
            covering = [
                seg
                for item in matched
                if item["reference_database_accession"] == accession
                for seg in item["aligned_regions"]
                if seg["entity_beg_seq_id"]
                <= position
                < seg["entity_beg_seq_id"] + seg["length"]
            ]
            if not covering:
                continue
            values = {
                seg["ref_beg_seq_id"] + position - seg["entity_beg_seq_id"]
                for seg in covering
            }
            if len(values) != 1:
                return base | {"mapping_reason": "conflicting_alignment"}
            value = values.pop()
            if sequence_entry and accession == sequence_entry["primary_accession"]:
                sequence = sequence_entry["sequence"]
                if value > len(sequence):
                    return base | {"mapping_reason": "uniprot_out_of_bounds"}
                base["uniprot_residue"] = sequence[value - 1]
            return base | {
                "uniprot_accession": accession,
                "uniprot_position": value,
                "mapping_reason": "mapped",
                "alignment_segments": covering,
            }
        return base | {"mapping_reason": "outside_alignment"}
    except RCSBError:
        return base | {"mapping_reason": "invalid_alignment"}


class Adapter:
    def __init__(self, accession, client, sleep, now, monotonic, uniprot_context):
        self.query = accession
        self.client, self.sleep, self.now, self.monotonic = (
            client,
            sleep,
            now,
            monotonic,
        )
        self.audit, self.cache = Audit(), {}
        self.last_status = None
        self.last_time = timestamp(now)
        self.accepted = [accession]
        self.entry = None
        if uniprot_context and (
            uniprot_context.get("primary_accession") == accession
            or accession in uniprot_context.get("secondary_accessions", [])
        ):
            self.entry = uniprot_context
            primary = uniprot_context["primary_accession"]
            if primary not in self.accepted:
                self.accepted.append(primary)
        self.preferred = self.entry["primary_accession"] if self.entry else accession
        self.context = {
            "query_accessions": self.accepted,
            "preferred_accession": self.preferred,
            "uniprot_resolution": {"query": accession, "primary": self.preferred}
            if self.entry
            else None,
            "coverage_limits": LIMITS.copy(),
            "entries": {},
            "target_entities": {},
            "ligand_entities": {},
        }
        self.selection = {"outcome_reason": None, "targets": [], "instances": []}
        self.warnings = []
        self.candidates = set()
        self.targets = {}
        self.records = {}

    def warn(self, code, message, role=None):
        warning = {"code": code, "message": message, "role": role}
        if warning not in self.warnings:
            self.warnings.append(warning)

    def request(self, kind, parts=(), params=None):
        url = (
            SEARCH
            if kind == "search"
            else CORE + "/".join((kind, *(quote(p, safe="") for p in parts)))
        )
        key = (url, json.dumps(params, sort_keys=True))
        if key in self.cache:
            return self.cache[key]
        roles = {
            "search": "search",
            "entry": "entry",
            "polymer_entity": "polymer_entity",
            "nonpolymer_entity": "nonpolymer_entity",
            "nonpolymer_entity_instance": "nonpolymer_instance",
            "chemcomp": "component",
        }
        response = request_json(
            url,
            params=params or {},
            role="rcsb_" + roles[kind],
            chemical_id=parts[0] if kind == "chemcomp" else None,
            audit=self.audit,
            client=self.client,
            sleep=self.sleep,
            now=self.now,
            monotonic=self.monotonic,
        )
        self.last_status, self.last_time = response.status, response.retrieved_at
        if response.error:
            raise RCSBError(response.error.message, response.error)
        value = (response.data, self.audit.request_count, response.status)
        self.cache[key] = value
        return value

    def search(self):
        result = set()
        for accession in self.accepted:
            start, total, seen = 0, None, set()
            while True:
                query = {
                    "query": {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                            "operator": "exact_match",
                            "value": accession,
                        },
                    },
                    "return_type": "polymer_entity",
                    "request_options": {
                        "results_content_type": ["experimental"],
                        "paginate": {"start": start, "rows": 100},
                        "sort": [{"sort_by": "rcsb_id", "direction": "asc"}],
                    },
                }
                data, _, status = self.request(
                    "search", params={"json": json.dumps(query, separators=(",", ":"))}
                )
                if status == 204:
                    demand(start == 0, "Search became empty during pagination.")
                    break
                demand(
                    isinstance(data, dict)
                    and data.get("result_type") == "polymer_entity"
                    and type(data.get("total_count")) is int
                    and data["total_count"] >= 0,
                    "Invalid RCSB search result.",
                )
                demand(
                    total is None or total == data["total_count"],
                    "Search total changed.",
                )
                total = data["total_count"]
                rows = array(data, "result_set")
                demand(
                    len(rows) <= 100
                    and start + len(rows) <= total
                    and (start + len(rows) == total or len(rows) == 100),
                    "Incomplete search page.",
                )
                for row in rows:
                    identifier = (
                        row.get("identifier") if isinstance(row, dict) else None
                    )
                    demand(
                        string(identifier) and "_" in identifier,
                        "Invalid entity search ID.",
                    )
                    entry, entity = identifier.rsplit("_", 1)
                    demand(
                        bool(re.fullmatch(r"[A-Za-z0-9]+", entry))
                        and bool(re.fullmatch(r"[1-9][0-9]*", entity)),
                        "Invalid entry/entity ID.",
                    )
                    demand(
                        identifier not in seen, "Duplicate entity in search pagination."
                    )
                    seen.add(identifier)
                start += len(rows)
                LOG.info("RCSB %s: %s/%s protein entities", accession, start, total)
                if start == total:
                    break
            result.update(seen)
        self.candidates = {i.rsplit("_", 1)[0] for i in result}
        return sorted(result)

    def object(self, kind, entry, entity=None):
        data, rid, _ = self.request(
            kind, (entry,) if entity is None else (entry, entity)
        )
        demand(isinstance(data, dict), f"Invalid {kind} object.")
        suffix = {
            "entry": "entry",
            "polymer_entity": "polymer_entity",
            "nonpolymer_entity": "nonpolymer_entity",
            "nonpolymer_entity_instance": "nonpolymer_entity_instance",
        }[kind]
        container = data.get("rcsb_" + suffix + "_container_identifiers")
        demand(
            isinstance(container, dict) and container.get("entry_id") == entry,
            f"{kind} parent mismatch.",
        )
        expected = (
            entry
            if entity is None
            else entry + ("." if kind.endswith("instance") else "_") + entity
        )
        demand(data.get("rcsb_id") == expected, f"{kind} identity mismatch.")
        if entity is not None:
            key = "asym_id" if kind.endswith("instance") else "entity_id"
            demand(container.get(key) == entity, f"{kind} child mismatch.")
        return data, container, rid

    def entry_metadata(self, entry):
        data, container, rid = self.object("entry", entry)
        info = data.get("rcsb_entry_info")
        demand(
            isinstance(info, dict)
            and string(info.get("structure_determination_methodology")),
            "Missing experimental methodology.",
        )
        if info["structure_determination_methodology"] != "experimental":
            return None
        entity_ids = identifiers(container, "non_polymer_entity_ids")
        demand(
            all(re.fullmatch(r"[1-9][0-9]*", i) for i in entity_ids),
            "Invalid nonpolymer entity IDs.",
        )
        meta = {"request_ids": [rid]}
        for key in ("struct", "exptl", "rcsb_primary_citation", "rcsb_accession_info"):
            if key in data:
                value = data[key]
                if isinstance(value, list if key == "exptl" else dict):
                    meta[key] = value
                else:
                    self.warn(
                        "invalid_optional_metadata",
                        f"{entry}: invalid {key}.",
                        "rcsb_entry",
                    )
        resolution = info.get("resolution_combined")
        if resolution is not None:
            if isinstance(resolution, list) and all(
                type(x) in (int, float) and math.isfinite(x) and x > 0
                for x in resolution
            ):
                meta["resolution_combined"] = resolution
            else:
                self.warn(
                    "invalid_optional_metadata",
                    f"{entry}: invalid resolution.",
                    "rcsb_entry",
                )
        self.context["entries"][entry] = meta
        return entity_ids

    def select_targets(self, identifiers_):
        entries = {}
        for identifier in identifiers_:
            entry, entity = identifier.rsplit("_", 1)
            data, container, rid = self.object("polymer_entity", entry, entity)
            chains = identifiers(container, "asym_ids", empty=False)
            refs = array(container, "reference_sequence_identifiers")
            demand(
                all(
                    isinstance(r, dict)
                    and string(r.get("database_name"))
                    and string(r.get("database_accession"))
                    for r in refs
                ),
                "Invalid reference identifiers.",
            )
            poly = data.get("entity_poly")
            demand(
                isinstance(poly, dict)
                and positive(poly.get("rcsb_sample_sequence_length")),
                "Invalid entity length.",
            )
            matching = [
                r
                for r in refs
                if r["database_name"] == "UniProt"
                and r["database_accession"] in self.accepted
            ]
            reason, request_ids = None, [rid]
            if not matching:
                reason = "accession_mismatch"
            else:
                if entry not in entries:
                    entries[entry] = self.entry_metadata(entry)
                if entries[entry] is None:
                    reason = "unsupported_methodology"
                else:
                    request_ids += self.context["entries"][entry]["request_ids"]
            self.selection["targets"].append(
                {
                    "entry_id": entry,
                    "entity_id": entity,
                    "included": reason is None,
                    "reason": reason,
                    "request_ids": sorted(request_ids),
                }
            )
            if reason:
                continue
            self.context["target_entities"][identifier] = {
                "asym_ids": chains,
                "references": refs,
                "sequence_length": poly["rcsb_sample_sequence_length"],
                "rcsb_polymer_entity": data.get("rcsb_polymer_entity", {}),
                "rcsb_polymer_entity_align": data.get("rcsb_polymer_entity_align", []),
                "request_ids": sorted(request_ids),
            }
            for chain in chains:
                demand(
                    (entry, chain) not in self.targets,
                    "Duplicate target chain assignment.",
                )
                self.targets[entry, chain] = (entity, data, rid, matching)
        return entries

    def contacts(self, entry, entity, ligand, comp, data, request_ids):
        rows = array(data, "rcsb_target_neighbors")
        groups = {}
        for row in rows:
            demand(
                isinstance(row, dict)
                and string(row.get("target_entity_id"))
                and string(row.get("target_asym_id")),
                "Invalid neighbor identity.",
            )
            chain = row["target_asym_id"]
            target = self.targets.get((entry, chain))
            if target is None:
                continue
            target_entity, target_data, target_rid, matching = target
            demand(
                row["target_entity_id"] == target_entity,
                "Neighbor target identity mismatch.",
            )
            demand(
                positive(row.get("target_seq_id"))
                and row["target_seq_id"]
                <= target_data["entity_poly"]["rcsb_sample_sequence_length"]
                and string(row.get("target_comp_id")),
                "Invalid target residue.",
            )
            distance = row.get("distance")
            demand(
                type(distance) in (int, float)
                and math.isfinite(distance)
                and distance >= 0,
                "Invalid contact distance.",
            )
            for key in ("target_model_id", "target_auth_seq_id"):
                if key in row:
                    demand(
                        type(row[key]) is int
                        and (key != "target_model_id" or row[key] > 0),
                        "Invalid neighbor numbering.",
                    )
            for key in ("atom_id", "target_atom_id", "alt_id", "target_is_bound"):
                demand(key not in row or isinstance(row[key], str), f"Invalid {key}.")
            demand(
                "comp_id" not in row or row["comp_id"] == comp,
                "Neighbor ligand mismatch.",
            )
            key = (chain, row.get("target_model_id"))
            groups.setdefault(key, []).append(row)
        sites, observations, accessions = [], [], set()
        container = data["rcsb_nonpolymer_entity_instance_container_identifiers"]
        for chain, model in sorted(
            groups, key=lambda k: (k[0], -1 if k[1] is None else k[1])
        ):
            target_entity, target_data, target_rid, matching = self.targets[
                entry, chain
            ]
            accessions.update(r["database_accession"] for r in matching)
            rows_ = groups[chain, model]
            residues = {}
            for row in rows_:
                pos = row["target_seq_id"]
                author = row.get("target_auth_seq_id")
                if pos in residues:
                    res = residues[pos]
                    demand(
                        res["comp_id"] == row["target_comp_id"]
                        and (
                            author is None
                            or res["auth_seq_id"] is None
                            or res["auth_seq_id"] == author
                        ),
                        "Conflicting residue identity.",
                    )
                    res["auth_seq_id"] = (
                        author if author is not None else res["auth_seq_id"]
                    )
                    res["distance"] = min(res["distance"], row["distance"])
                    continue
                mapped = mapping(
                    target_data, pos, self.accepted, self.preferred, self.entry
                )
                if mapped["uniprot_residue"] is not None:
                    observed = AMINO.get(row["target_comp_id"])
                    mapped["sequence_check"] = (
                        "unknown_residue"
                        if observed is None
                        else "match"
                        if observed == mapped["uniprot_residue"]
                        else "mismatch"
                    )
                if (
                    mapped["mapping_reason"]
                    in (
                        "invalid_alignment",
                        "conflicting_alignment",
                        "uniprot_out_of_bounds",
                    )
                    or mapped["sequence_check"] == "mismatch"
                ):
                    self.warn(
                        "residue_mapping",
                        f"{entry}.{chain} residue {pos}: {mapped['mapping_reason']}, {mapped['sequence_check']}.",
                        "rcsb_polymer_entity",
                    )
                residues[pos] = dict(
                    label_seq_id=pos,
                    auth_seq_id=author,
                    comp_id=row["target_comp_id"],
                    distance=row["distance"],
                    **mapped,
                )
            ids = sorted(set(request_ids + [target_rid]))
            identity = {
                "entry_id": entry,
                "ligand_entity_id": entity,
                "ligand_asym_id": ligand,
                "target_entity_id": target_entity,
                "target_asym_id": chain,
                "target_model_id": model,
            }
            sites.append(
                identity
                | {
                    "ligand_auth_asym_id": container.get("auth_asym_id"),
                    "ligand_auth_seq_id": container.get("auth_seq_id"),
                    "target_auth_asym_id": None,
                    "residues": [residues[p] for p in sorted(residues)],
                    "request_ids": ids,
                }
            )
            observations.append(identity | {"neighbor_rows": rows_, "request_ids": ids})
        return sites, observations, accessions, bool(rows)

    def ligands(self, entries):
        for entry in sorted(entries):
            if entries[entry] is None:
                continue
            entry_request = self.context["entries"][entry]["request_ids"]
            for entity in sorted(entries[entry]):
                data, container, rid = self.object("nonpolymer_entity", entry, entity)
                chains = identifiers(container, "asym_ids", empty=False)
                comp = container.get("nonpolymer_comp_id")
                demand(string(comp), "Missing ligand CCD identity.")
                annotation = data.get("pdbx_entity_nonpoly", {})
                demand(
                    isinstance(annotation, dict)
                    and annotation.get("comp_id", comp) == comp,
                    "Ligand chemical identity mismatch.",
                )
                self.context["ligand_entities"][f"{entry}_{entity}"] = {
                    "ccd_id": comp,
                    "name": annotation.get("name")
                    if string(annotation.get("name"))
                    else None,
                    "request_ids": [rid],
                }
                for chain in sorted(chains):
                    instance, ci, iid = self.object(
                        "nonpolymer_entity_instance", entry, chain
                    )
                    demand(
                        ci.get("entity_id") == entity and ci.get("comp_id") == comp,
                        "Ligand instance parent mismatch.",
                    )
                    for key in ("auth_asym_id", "auth_seq_id"):
                        demand(
                            key not in ci or isinstance(ci[key], str),
                            "Invalid ligand author numbering.",
                        )
                    ids = sorted(set(entry_request + [rid, iid]))
                    sites, obs, acc, has_neighbors = (
                        ([], [], set(), False)
                        if comp in ("HOH", "DOD")
                        else self.contacts(entry, entity, chain, comp, instance, ids)
                    )
                    reason = (
                        "water"
                        if comp in ("HOH", "DOD")
                        else None
                        if sites
                        else "no_target_contact"
                        if has_neighbors
                        else "no_reported_contacts"
                    )
                    self.selection["instances"].append(
                        {
                            "entry_id": entry,
                            "entity_id": entity,
                            "asym_id": chain,
                            "ccd_id": comp,
                            "included": reason is None,
                            "reason": reason,
                            "request_ids": ids,
                        }
                    )
                    if reason:
                        continue
                    if comp not in self.records:
                        self.records[comp] = Measurement(
                            compound_id="CCD:" + comp,
                            compound_name=annotation.get("name")
                            if string(annotation.get("name"))
                            else comp,
                            smiles=None,
                            inchikey=None,
                            target_description=None,
                            evidence_type="structure_associated_ligand",
                            measurement_type=None,
                            measurement_value=None,
                            measurement_unit=None,
                            pmid=None,
                            doi=None,
                            source_url="https://www.rcsb.org/ligand/"
                            + quote(comp, safe=""),
                            target_match_status="exact_uniprot_reference",
                            source_values={"observations": []},
                        )
                    record = self.records[comp]
                    record.binding_sites.extend(sites)
                    record.source_values["observations"].extend(obs)
                    record.target_accessions = sorted(
                        set(record.target_accessions) | acc
                    )
            LOG.info("RCSB %s: processed entry %s", self.query, entry)

    def enrich(self):
        from ligandscout.ccd import descriptors

        for comp, record in sorted(self.records.items()):

            def site_key(item):
                return (
                    item["entry_id"],
                    item["ligand_asym_id"],
                    item["target_asym_id"],
                    -1 if item["target_model_id"] is None else item["target_model_id"],
                )

            record.binding_sites.sort(key=site_key)
            record.source_values["observations"].sort(key=site_key)
            prov = {
                "provider": "RCSB CCD",
                "status": "unavailable",
                "ccd_id": comp,
                "smiles_field": None,
                "request_ids": [],
            }
            try:
                data, rid, _ = self.request("chemcomp", (comp,))
                prov["request_ids"] = [rid]
                demand(
                    isinstance(data, dict)
                    and data.get("rcsb_id") == comp
                    and isinstance(data.get("chem_comp"), dict)
                    and data["chem_comp"].get("id") == comp,
                    "CCD component identity mismatch.",
                )
                field, smiles, key = descriptors(data)
                record.smiles, record.inchikey = smiles, key
                if string(data["chem_comp"].get("name")):
                    record.compound_name = data["chem_comp"]["name"]
                prov.update(
                    status="available" if smiles is not None else "no_descriptor",
                    smiles_field=field,
                )
            except (RCSBError, ValueError) as exc:
                prov["request_ids"] = [self.audit.request_count]
                self.warn("structure_unavailable", f"{comp}: {exc}", "rcsb_component")
            if prov["status"] == "no_descriptor":
                self.warn(
                    "structure_no_descriptor",
                    f"{comp}: SMILES unavailable.",
                    "rcsb_component",
                )
            record.structure_provenance = prov


def fetch_rcsb(accession, *, client, sleep, now, monotonic, uniprot_context=None):
    accession = normalize_accession(accession)
    adapter = Adapter(accession, client, sleep, now, monotonic, uniprot_context)
    error = None
    LOG.info("Searching RCSB PDB for %s", accession)
    try:
        entities = adapter.search()
        entries = adapter.select_targets(entities)
        adapter.ligands(entries)
        adapter.enrich()
        records = [adapter.records[k] for k in sorted(adapter.records)]
        adapter.selection["outcome_reason"] = (
            None
            if records
            else "no_exact_structures"
            if not entities
            else "no_eligible_entities"
            if not adapter.targets
            else "no_reported_target_ligands"
        )
    except RCSBError as exc:
        error = exc.error or SourceError("schema", str(exc), adapter.last_status)
        records = []
        adapter.selection["outcome_reason"] = "incomplete"
    result = SourceResult(
        schema_version=5,
        accession=accession,
        source="RCSB PDB",
        retrieved_at=adapter.last_time,
        request={
            "search_endpoint": SEARCH,
            "query_accessions": adapter.accepted,
            "page_size": 100,
            "return_type": "polymer_entity",
            "results_content_type": ["experimental"],
        },
        source_version=None,
        status="failed" if error else "success" if records else "no_matches",
        records=records,
        error=error,
        policy="rcsb-structural-contacts-v1",
        warnings=adapter.warnings,
        counts=counts(
            records, adapter.selection, adapter.candidates, failed=error is not None
        ),
        selection_audit=adapter.selection,
        source_context=adapter.context,
    )
    validate_result(result)
    LOG.info(
        "RCSB %s: %s; %s ligands, %s sites",
        accession,
        result.status,
        result.counts["unique_ligands"],
        result.counts["contact_sites"],
    )
    return SearchOutcome(result, None, adapter.audit.exchanges)
