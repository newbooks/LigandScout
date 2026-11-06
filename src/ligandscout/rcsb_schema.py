"""Structural ligand count and saved-result invariants."""

import math

from ligandscout.models import ANNOTATION_COUNTS, SOURCE_POLICIES
from ligandscout.uniprot_schema import ACTIVITY_COUNTS

STRUCTURAL_COUNTS = (
    "candidate_entries",
    "accepted_target_entities",
    "examined_ligand_instances",
    "included_ligand_instances",
    "excluded_water_instances",
    "instances_without_reported_contacts",
    "instances_without_target_contacts",
    "contact_sites",
    "contact_residues",
    "mapped_contact_residues",
)


def counts(records, audit, candidates, *, failed=False):
    result = dict.fromkeys((*ACTIVITY_COUNTS, *ANNOTATION_COUNTS))
    result.update(complete=not failed, record_kind="structural_ligand")
    keys = (
        *STRUCTURAL_COUNTS,
        "unique_ligands",
        "smiles_available",
        "smiles_unavailable",
    )
    if failed:
        return result | dict.fromkeys(keys)
    sites = [s for r in records for s in r.binding_sites]
    residues = [v for s in sites for v in s["residues"]]
    instances = audit["instances"]
    result.update(
        candidate_entries=len(candidates),
        accepted_target_entities=sum(t["included"] for t in audit["targets"]),
        examined_ligand_instances=len(instances),
        included_ligand_instances=sum(i["included"] for i in instances),
        excluded_water_instances=sum(i["reason"] == "water" for i in instances),
        instances_without_reported_contacts=sum(
            i["reason"] == "no_reported_contacts" for i in instances
        ),
        instances_without_target_contacts=sum(
            i["reason"] == "no_target_contact" for i in instances
        ),
        contact_sites=len(sites),
        contact_residues=len(residues),
        mapped_contact_residues=sum(
            r["uniprot_position"] is not None for r in residues
        ),
        unique_ligands=len(records),
        smiles_available=sum(r.smiles is not None for r in records),
        smiles_unavailable=sum(r.smiles is None for r in records),
    )
    return result


def require(condition, message):
    if not condition:
        raise ValueError(message)


def _validate_result(result):
    require(
        result.schema_version == 5 and result.policy == SOURCE_POLICIES["RCSB PDB"],
        "Invalid RCSB schema/policy.",
    )
    require(
        result.status in ("success", "no_matches", "failed"), "Invalid RCSB status."
    )
    require(
        (result.status == "failed") == (result.error is not None), "Invalid RCSB error."
    )
    require(
        (result.status == "success") == bool(result.records),
        "Invalid RCSB records/status.",
    )
    require(
        set(result.selection_audit) == {"outcome_reason", "targets", "instances"},
        "Invalid RCSB audit.",
    )
    context = result.source_context
    require(
        set(context)
        == {
            "query_accessions",
            "preferred_accession",
            "uniprot_resolution",
            "coverage_limits",
            "entries",
            "target_entities",
            "ligand_entities",
        },
        "Invalid RCSB context.",
    )
    for warning in result.warnings:
        require(
            set(warning) == {"code", "message", "role"}
            and isinstance(warning["code"], str)
            and isinstance(warning["message"], str)
            and (warning["role"] is None or isinstance(warning["role"], str)),
            "Invalid warning.",
        )
    require(
        result.counts.get("record_kind") == "structural_ligand", "Invalid record kind."
    )
    expected_keys = {
        "complete",
        "record_kind",
        "unique_ligands",
        *ACTIVITY_COUNTS,
        *ANNOTATION_COUNTS,
        *STRUCTURAL_COUNTS,
    }
    require(set(result.counts) == expected_keys, "Invalid structural count keys.")
    for k, v in result.counts.items():
        if k not in ("complete", "record_kind"):
            require(v is None or type(v) is int and v >= 0, "Invalid count value.")
    for kind in ("targets", "instances"):
        require(isinstance(result.selection_audit[kind], list), "Invalid audit rows.")
        keys = {"entry_id", "entity_id", "included", "reason", "request_ids"}
        if kind == "instances":
            keys |= {"asym_id", "ccd_id"}
        for row in result.selection_audit[kind]:
            require(isinstance(row, dict) and set(row) == keys, "Invalid audit row.")
            require(
                type(row["included"]) is bool
                and (
                    row["reason"] is None
                    if row["included"]
                    else row["reason"]
                    in {
                        "accession_mismatch",
                        "unsupported_methodology",
                        "water",
                        "no_reported_contacts",
                        "no_target_contact",
                    }
                ),
                "Invalid exclusion reason.",
            )
    for field in ("entries", "target_entities", "ligand_entities"):
        require(isinstance(context[field], dict), "Invalid context object.")
    require(
        isinstance(context["query_accessions"], list)
        and context["query_accessions"]
        and all(isinstance(x, str) for x in context["query_accessions"])
        and result.accession == context["query_accessions"][0]
        and context["preferred_accession"] in context["query_accessions"],
        "Invalid query accessions.",
    )
    require(type(result.counts["complete"]) is bool, "Invalid completeness.")
    for kind in ("targets", "instances"):
        audit_ids = []
        for row in result.selection_audit[kind]:
            require(
                all(
                    isinstance(row[k], str) and row[k]
                    for k in ("entry_id", "entity_id")
                ),
                "Invalid audit identity.",
            )
            require(
                isinstance(row["request_ids"], list)
                and row["request_ids"]
                and all(type(x) is int and x > 0 for x in row["request_ids"]),
                "Invalid audit requests.",
            )
            audit_ids.append(
                (row["entry_id"], row["entity_id"])
                if kind == "targets"
                else (row["entry_id"], row["asym_id"])
            )
        require(len(set(audit_ids)) == len(audit_ids), "Duplicate audit identity.")
    ids = set()
    for record in result.records:
        require(
            record.compound_id.startswith("CCD:") and record.compound_id not in ids,
            "Invalid/duplicate CCD group.",
        )
        ids.add(record.compound_id)
        require(
            record.evidence_type == "structure_associated_ligand"
            and record.target_match_status == "exact_uniprot_reference",
            "Invalid structural evidence.",
        )
        require(
            all(
                getattr(record, k) is None
                for k in (
                    "measurement_type",
                    "measurement_value",
                    "measurement_unit",
                    "pmid",
                    "doi",
                    "assay_id",
                    "target_id",
                )
            ),
            "Structural record contains affinity data.",
        )
        require(
            record.ligand_roles == [] and record.annotation_evidence == [],
            "Invalid structural annotation evidence.",
        )
        require(
            record.binding_sites
            and len(record.binding_sites) == len(record.source_values["observations"]),
            "Invalid sites/observations.",
        )
        seen = set()
        for site, obs in zip(
            record.binding_sites, record.source_values["observations"], strict=True
        ):
            site_keys = {
                "entry_id",
                "ligand_entity_id",
                "ligand_asym_id",
                "ligand_auth_asym_id",
                "ligand_auth_seq_id",
                "target_entity_id",
                "target_asym_id",
                "target_auth_asym_id",
                "target_model_id",
                "residues",
                "request_ids",
            }
            obs_keys = {
                "entry_id",
                "ligand_entity_id",
                "ligand_asym_id",
                "target_entity_id",
                "target_asym_id",
                "target_model_id",
                "neighbor_rows",
                "request_ids",
            }
            require(
                isinstance(site, dict)
                and set(site) == site_keys
                and isinstance(obs, dict)
                and set(obs) == obs_keys,
                "Invalid site/observation fields.",
            )
            require(
                site["target_auth_asym_id"] is None
                and (
                    site["target_model_id"] is None
                    or type(site["target_model_id"]) is int
                    and site["target_model_id"] > 0
                ),
                "Invalid target model/author chain.",
            )
            require(
                all(
                    site[k] is None or isinstance(site[k], str)
                    for k in ("ligand_auth_asym_id", "ligand_auth_seq_id")
                ),
                "Invalid ligand numbering.",
            )
            require(
                site["request_ids"] == obs["request_ids"]
                and isinstance(site["request_ids"], list)
                and site["request_ids"]
                and all(type(i) is int and i > 0 for i in site["request_ids"]),
                "Invalid site requests.",
            )
            key = tuple(
                site[k]
                for k in (
                    "entry_id",
                    "ligand_asym_id",
                    "target_asym_id",
                    "target_model_id",
                )
            )
            require(
                key not in seen
                and all(
                    obs[k] == site[k]
                    for k in (
                        "entry_id",
                        "ligand_entity_id",
                        "ligand_asym_id",
                        "target_entity_id",
                        "target_asym_id",
                        "target_model_id",
                    )
                ),
                "Invalid structural site identity.",
            )
            seen.add(key)
            entry, entity = site["entry_id"], site["target_entity_id"]
            require(
                entry in context["entries"]
                and f"{entry}_{entity}" in context["target_entities"],
                "Missing structural context.",
            )
            target = context["target_entities"][f"{entry}_{entity}"]
            require(
                site["target_asym_id"] in target["asym_ids"], "Unknown target chain."
            )
            require(obs["neighbor_rows"] and site["residues"], "Empty site.")
            require(
                isinstance(site["residues"], list)
                and isinstance(obs["neighbor_rows"], list),
                "Invalid contact collections.",
            )
            ligand_key = f"{entry}_{site['ligand_entity_id']}"
            require(
                ligand_key in context["ligand_entities"]
                and context["ligand_entities"][ligand_key]["ccd_id"]
                == record.compound_id[4:],
                "Invalid ligand parent context.",
            )
            audit_instance = [
                i
                for i in result.selection_audit["instances"]
                if i["entry_id"] == entry
                and i["asym_id"] == site["ligand_asym_id"]
                and i["included"]
            ]
            require(
                len(audit_instance) == 1
                and audit_instance[0]["entity_id"] == site["ligand_entity_id"]
                and audit_instance[0]["ccd_id"] == record.compound_id[4:],
                "Site missing included instance.",
            )
            positions = set()
            for residue in site["residues"]:
                residue_keys = {
                    "label_seq_id",
                    "auth_seq_id",
                    "comp_id",
                    "distance",
                    "uniprot_accession",
                    "uniprot_position",
                    "mapping_reason",
                    "uniprot_residue",
                    "sequence_check",
                    "alignment_segments",
                }
                require(
                    isinstance(residue, dict) and set(residue) == residue_keys,
                    "Invalid residue fields.",
                )
                require(
                    isinstance(residue["comp_id"], str)
                    and residue["comp_id"]
                    and type(residue["distance"]) in (int, float)
                    and math.isfinite(residue["distance"])
                    and residue["distance"] >= 0,
                    "Invalid residue identity/distance.",
                )
                require(
                    residue["auth_seq_id"] is None
                    or type(residue["auth_seq_id"]) is int,
                    "Invalid author position.",
                )
                require(
                    residue["mapping_reason"]
                    in {
                        "mapped",
                        "no_alignment",
                        "outside_alignment",
                        "conflicting_alignment",
                        "invalid_alignment",
                        "uniprot_out_of_bounds",
                    }
                    and residue["sequence_check"]
                    in {"match", "mismatch", "unknown_residue", "not_checked"},
                    "Invalid mapping reason.",
                )
                pos = residue["label_seq_id"]
                require(
                    type(pos) is int
                    and 1 <= pos <= target["sequence_length"]
                    and pos not in positions,
                    "Invalid residue position.",
                )
                positions.add(pos)
                source_rows = [
                    row
                    for row in obs["neighbor_rows"]
                    if row.get("target_seq_id") == pos
                ]
                require(
                    source_rows
                    and all(
                        row["target_entity_id"] == entity
                        and row["target_asym_id"] == site["target_asym_id"]
                        and row.get("target_model_id") == site["target_model_id"]
                        and row["target_comp_id"] == residue["comp_id"]
                        for row in source_rows
                    ),
                    "Residue disagrees with neighbor source.",
                )
                require(
                    residue["distance"] == min(row["distance"] for row in source_rows),
                    "Residue distance disagrees with source.",
                )
                mapped = residue["uniprot_position"]
                require(
                    (mapped is not None) == (residue["mapping_reason"] == "mapped"),
                    "Invalid mapped residue.",
                )
                require(
                    mapped is None
                    or type(mapped) is int
                    and mapped > 0
                    and residue["uniprot_accession"] in context["query_accessions"],
                    "Invalid UniProt position.",
                )
            require(
                positions == {row["target_seq_id"] for row in obs["neighbor_rows"]},
                "Missing source residues.",
            )
        prov = record.structure_provenance
        require(
            set(prov) == {"provider", "status", "ccd_id", "smiles_field", "request_ids"}
            and prov["provider"] == "RCSB CCD"
            and prov["ccd_id"] == record.compound_id[4:],
            "Invalid CCD provenance.",
        )
        require(
            prov["status"] in ("available", "no_descriptor", "unavailable")
            and (prov["status"] == "available") == (record.smiles is not None),
            "Invalid CCD outcome.",
        )
    if result.status == "failed":
        require(
            result.counts["complete"] is False
            and all(
                v is None
                for k, v in result.counts.items()
                if k not in ("complete", "record_kind")
            ),
            "Failed source has totals.",
        )
    else:
        # Candidate IDs are deliberately recorded separately from fetched entries.
        candidate_count = result.counts["candidate_entries"]
        require(
            type(candidate_count) is int and candidate_count >= 0,
            "Invalid candidate count.",
        )
        expected = counts(
            result.records, result.selection_audit, range(candidate_count)
        )
        require(result.counts == expected, "Structural counts disagree with evidence.")
        c = result.counts
        require(
            c["examined_ligand_instances"]
            == sum(
                c[k]
                for k in (
                    "included_ligand_instances",
                    "excluded_water_instances",
                    "instances_without_reported_contacts",
                    "instances_without_target_contacts",
                )
            ),
            "Invalid instance partition.",
        )
        require(
            c["unique_ligands"] <= c["included_ligand_instances"]
            and c["mapped_contact_residues"] <= c["contact_residues"],
            "Invalid count bounds.",
        )


def validate_result(result):
    """Reject malformed saved structural data with a consistent validation error."""
    try:
        _validate_result(result)
    except (KeyError, TypeError, AttributeError, OverflowError) as exc:
        raise ValueError("Malformed structural result.") from exc
