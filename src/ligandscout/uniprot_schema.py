"""Schema-v4 annotation counts and saved-result invariants."""

import hashlib
import re

from ligandscout.accessions import normalize_accession
from ligandscout.models import ANNOTATION_COUNTS, SOURCE_POLICIES

ACTIVITY_COUNTS = (
    "candidate_targets",
    "eligible_targets",
    "candidate_activities",
    "included_activities",
    "excluded_activities",
    "unverified_activities",
    "candidate_assays",
    "eligible_assays",
)


class UniProtSchemaError(ValueError):
    """UniProt annotation data is inconsistent or invalid."""


def annotation_counts(records, context, *, failed=False):
    base = dict.fromkeys(ACTIVITY_COUNTS)
    return (
        base
        | {"complete": not failed, "record_kind": "ligand_annotation"}
        | (
            dict.fromkeys((*ANNOTATION_COUNTS, "unique_ligands"))
            if failed
            else {
                "unique_ligands": len(records),
                "annotation_records": len(records),
                "identified_chemical_ids": sum(
                    r.compound_id.startswith("CHEBI:") for r in records
                ),
                "binding_site_annotations": sum(len(r.binding_sites) for r in records),
                "generic_binding_sites": len(context.get("generic_sites", [])),
                "cofactor_annotations": sum(
                    len(r.source_values["cofactor_annotations"]) for r in records
                )
                + len(context.get("generic_cofactors", [])),
                "smiles_available": sum(bool(r.smiles) for r in records),
                "smiles_unavailable": sum(not r.smiles for r in records),
            }
        )
    )


def demand(condition, message):
    if not condition:
        raise ValueError(message)


def validate_evidence(items):
    for item in items:
        demand(
            isinstance(item, dict)
            and set(item)
            == {"annotation_path", "evidence_code", "source", "source_id", "category"},
            "Invalid evidence fields.",
        )
        demand(
            isinstance(item["annotation_path"], str)
            and item["annotation_path"].startswith("/"),
            "Invalid evidence path.",
        )
        demand(
            isinstance(item["evidence_code"], str)
            and re.fullmatch(r"ECO:[0-9]{7}", item["evidence_code"]),
            "Invalid evidence code.",
        )
        demand(
            all(
                item[k] is None or isinstance(item[k], str)
                for k in ("source", "source_id")
            ),
            "Invalid evidence source.",
        )
        expected = {
            "ECO:0000269": "experimental",
            "ECO:0000250": "similarity",
            "ECO:0000255": "sequence_analysis",
        }.get(item["evidence_code"], "other")
        demand(item["category"] == expected, "Invalid evidence category.")


def validate_sites(sites, entry):
    for site in sites:
        demand(
            set(site)
            == {
                "annotation_path",
                "feature_id",
                "location",
                "sequence_scope",
                "residue_sequence",
                "description",
                "ligand_part",
                "ligand_label",
                "processed_products",
                "evidence",
            },
            "Invalid site fields.",
        )
        demand(isinstance(site["annotation_path"], str), "Invalid site path.")
        demand(
            all(
                site[k] is None or isinstance(site[k], str)
                for k in (
                    "feature_id",
                    "residue_sequence",
                    "description",
                    "ligand_label",
                )
            ),
            "Invalid site text.",
        )
        location = site["location"]
        demand(
            isinstance(location, dict)
            and set(location) == {"start", "end", "sequence"},
            "Invalid location.",
        )
        for key in ("start", "end"):
            value = location[key]
            demand(
                isinstance(value, dict) and set(value) == {"value", "modifier"},
                "Invalid endpoint.",
            )
            demand(
                isinstance(value["modifier"], str) and bool(value["modifier"]),
                "Invalid coordinate modifier.",
            )
            demand(
                value["value"] is None
                or type(value["value"]) is int
                and value["value"] > 0,
                "Invalid coordinate.",
            )
            demand(
                value["modifier"] != "EXACT" or value["value"] is not None,
                "Missing exact coordinate.",
            )
        demand(
            location["sequence"] is None or isinstance(location["sequence"], str),
            "Invalid sequence ID.",
        )
        canonical = location["sequence"] in (None, entry["primary_accession"])
        demand(
            site["sequence_scope"] == ("canonical" if canonical else "alternate"),
            "Incorrect site sequence scope.",
        )
        exact = location["start"]["modifier"] == location["end"]["modifier"] == "EXACT"
        expected = None
        if exact:
            start, end = location["start"]["value"], location["end"]["value"]
            demand(start <= end, "Reversed site.")
            if canonical:
                demand(end <= entry["sequence_length"], "Out-of-bounds site.")
                expected = entry["sequence"][start - 1 : end]
        demand(site["residue_sequence"] == expected, "Incorrect site residues.")
        demand(
            isinstance(site["processed_products"], list)
            and isinstance(site["evidence"], list),
            "Invalid site arrays.",
        )
        part = site["ligand_part"]
        demand(
            part is None
            or isinstance(part, dict)
            and set(part) == {"name", "id"}
            and all(x is None or isinstance(x, str) for x in part.values()),
            "Invalid ligand part.",
        )
        for product in site["processed_products"]:
            demand(
                isinstance(product, dict)
                and set(product) == {"annotation_path", "name", "start", "end"},
                "Invalid processed product fields.",
            )
            demand(
                isinstance(product["annotation_path"], str)
                and (product["name"] is None or isinstance(product["name"], str)),
                "Invalid processed product name/path.",
            )
            demand(
                type(product["start"]) is int
                and type(product["end"]) is int
                and 1 <= product["start"] <= product["end"] <= entry["sequence_length"],
                "Invalid processed product coordinates.",
            )
            demand(
                canonical
                and exact
                and product["start"] <= start
                and product["end"] >= end,
                "Product does not contain the site.",
            )
        validate_evidence(site["evidence"])


def validate_result(result):
    demand(
        result.schema_version == 4 and result.policy == SOURCE_POLICIES["UniProt"],
        "Invalid UniProt schema/policy.",
    )
    demand(
        result.counts.get("record_kind") == "ligand_annotation", "Invalid record kind."
    )
    for warning in result.warnings:
        demand(
            set(warning) == {"code", "message", "role"}
            and isinstance(warning["code"], str)
            and isinstance(warning["message"], str)
            and (warning["role"] is None or isinstance(warning["role"], str)),
            "Invalid saved warning.",
        )
    if result.status == "failed":
        demand(
            result.source_context == {} and result.records == [],
            "Failed source has partial annotations.",
        )
    else:
        context = result.source_context
        demand(
            set(context)
            == {
                "entry",
                "generic_sites",
                "generic_cofactors",
                "context_comments",
                "processed_products",
            },
            "Invalid UniProt context.",
        )
        entry = context["entry"]
        demand(
            isinstance(entry, dict)
            and isinstance(entry.get("primary_accession"), str)
            and isinstance(entry.get("sequence"), str)
            and type(entry.get("sequence_length")) is int
            and len(entry["sequence"]) == entry["sequence_length"],
            "Invalid saved entry sequence.",
        )
        demand(
            all(
                isinstance(context[k], list)
                for k in (
                    "generic_sites",
                    "generic_cofactors",
                    "context_comments",
                    "processed_products",
                )
            ),
            "Invalid context collections.",
        )
        demand(
            entry.get("query_accession") == result.accession
            and isinstance(entry.get("secondary_accessions"), list),
            "Invalid saved query identity.",
        )
        demand(
            normalize_accession(entry["primary_accession"])
            == entry["primary_accession"]
            and (
                result.accession == entry["primary_accession"]
                or result.accession in entry["secondary_accessions"]
            ),
            "Invalid resolved accession.",
        )
        demand(
            entry.get("entry_type") == "UniProtKB reviewed (Swiss-Prot)"
            and re.fullmatch(r"[A-Z]+", entry["sequence"])
            and isinstance(entry.get("organism"), str)
            and type(entry.get("taxon_id")) is int
            and entry["taxon_id"] > 0,
            "Invalid saved reviewed entry.",
        )
        for comment in context["context_comments"]:
            demand(
                isinstance(comment, dict)
                and set(comment) == {"comment_type", "annotation_path", "data"}
                and isinstance(comment["comment_type"], str)
                and isinstance(comment["annotation_path"], str)
                and isinstance(comment["data"], dict),
                "Invalid saved context comment.",
            )
        validate_sites(context["generic_sites"], entry)
        ids = []
        for record in result.records:
            demand(record.compound_id not in ids, "Duplicate UniProt ligand group.")
            ids.append(record.compound_id)
            demand(
                record.evidence_type == "uniprot_annotation"
                and record.target_match_status == "entry_annotation"
                and record.target_id == entry["primary_accession"]
                and record.target_accessions == [record.target_id],
                "Invalid annotation target.",
            )
            demand(
                record.ligand_roles
                and record.ligand_roles
                == sorted(
                    set(record.ligand_roles),
                    key=lambda v: ("binding_site", "cofactor").index(v),
                ),
                "Invalid ligand roles.",
            )
            demand(
                all(
                    getattr(record, k) is None
                    for k in (
                        "measurement_type",
                        "measurement_value",
                        "measurement_unit",
                        "measurement_relation",
                        "measurement_upper_value",
                        "measurement_text",
                        "measurement_basis",
                        "pmid",
                        "doi",
                        "assay_id",
                        "assay_type",
                        "source_record_id",
                        "document_id",
                        "target_confidence",
                        "target_relationship",
                        "variant_accession",
                        "variant_mutation",
                        "variant_isoform",
                    )
                ),
                "Annotations contain fabricated assay metadata.",
            )
            demand(
                set(record.source_values)
                == {"aliases", "binding_annotations", "cofactor_annotations"}
                and all(isinstance(v, list) for v in record.source_values.values()),
                "Invalid annotation provenance.",
            )
            aliases = record.source_values["aliases"]
            demand(
                all(isinstance(name, str) and name.strip() for name in aliases),
                "Invalid ligand aliases.",
            )
            if record.compound_id.startswith("CHEBI:"):
                demand(
                    re.fullmatch(r"CHEBI:[1-9][0-9]*", record.compound_id),
                    "Invalid chemical ID.",
                )
            else:
                demand(bool(aliases), "Name-only ligand lacks a name.")
                digest = hashlib.sha256(
                    " ".join(aliases[0].split()).casefold().encode()
                ).hexdigest()
                demand(
                    record.compound_id
                    == f"UNIPROT:{entry['primary_accession']}:NAME:{digest}",
                    "Invalid local ligand ID.",
                )
            demand(
                record.source_url
                == f"https://www.uniprot.org/uniprotkb/{record.target_id}/entry",
                "Invalid annotation source URL.",
            )
            demand(
                len(record.source_values["binding_annotations"])
                == len(record.binding_sites),
                "Binding annotation/site counts differ.",
            )
            expected_roles = [
                role
                for role, values in (
                    ("binding_site", record.binding_sites),
                    ("cofactor", record.source_values["cofactor_annotations"]),
                )
                if values
            ]
            demand(
                record.ligand_roles == expected_roles,
                "Roles differ from source annotations.",
            )
            for annotation in (
                record.source_values["binding_annotations"]
                + record.source_values["cofactor_annotations"]
            ):
                demand(
                    isinstance(annotation, dict)
                    and isinstance(annotation.get("annotation_path"), str)
                    and isinstance(annotation.get("data"), dict),
                    "Invalid original annotation.",
                )
            validate_sites(record.binding_sites, entry)
            validate_evidence(record.annotation_evidence)
            structure = record.structure_provenance
            demand(
                set(structure)
                == {
                    "provider",
                    "status",
                    "requested_chemical_id",
                    "smiles_field",
                    "candidate_components",
                    "request_ids",
                },
                "Invalid structure provenance.",
            )
            demand(
                structure["provider"] == "RCSB CCD"
                and structure["status"]
                in (
                    "available",
                    "no_exact_mapping",
                    "no_descriptor",
                    "ambiguous",
                    "unavailable",
                    "not_applicable",
                ),
                "Invalid structure status.",
            )
            demand(
                isinstance(structure["candidate_components"], list)
                and isinstance(structure["request_ids"], list),
                "Invalid structure candidates/requests.",
            )
            demand(
                all(type(i) is int and i > 0 for i in structure["request_ids"]),
                "Invalid structure request IDs.",
            )
            chemical = (
                record.compound_id if record.compound_id.startswith("CHEBI:") else None
            )
            demand(
                structure["requested_chemical_id"] == chemical,
                "Structure chemical ID mismatch.",
            )
            demand(
                structure["status"] != "not_applicable" or chemical is None,
                "Applicable chemical lookup was skipped.",
            )
            for candidate in structure["candidate_components"]:
                demand(
                    isinstance(candidate, dict)
                    and set(candidate)
                    == {
                        "ccd_id",
                        "mapping",
                        "smiles_field",
                        "smiles",
                        "inchikey",
                        "request_id",
                        "reason",
                    },
                    "Invalid CCD candidate fields.",
                )
                demand(
                    isinstance(candidate["ccd_id"], str)
                    and re.fullmatch(r"[A-Za-z0-9]+", candidate["ccd_id"])
                    and type(candidate["request_id"]) is int
                    and candidate["request_id"] in structure["request_ids"],
                    "Invalid CCD candidate identity/request.",
                )
                demand(
                    all(
                        candidate[k] is None or isinstance(candidate[k], str)
                        for k in ("smiles_field", "smiles", "inchikey", "reason")
                    ),
                    "Invalid CCD descriptor text.",
                )
                demand(
                    candidate["mapping"] is None
                    or isinstance(candidate["mapping"], dict),
                    "Invalid CCD mapping.",
                )
            if structure["status"] == "available":
                demand(
                    record.smiles
                    and structure["smiles_field"] in ("SMILES", "SMILES_stereo"),
                    "Missing available structure.",
                )
                eligible = [
                    c for c in structure["candidate_components"] if c.get("mapping")
                ]
                demand(
                    eligible
                    and all(
                        (c.get("smiles"), c.get("inchikey"), c.get("smiles_field"))
                        == (record.smiles, record.inchikey, structure["smiles_field"])
                        and c["mapping"].get("resource_name") == "ChEBI"
                        and c["mapping"].get("resource_accession_code")
                        == record.compound_id
                        for c in eligible
                    ),
                    "Structure does not match exact chemical identifier.",
                )
            else:
                demand(
                    record.smiles is None and record.inchikey is None,
                    "Unresolved structure has a selected SMILES.",
                )
        demand(
            result.status == ("success" if result.records else "no_matches"),
            "Inconsistent annotation status.",
        )
    expected = annotation_counts(
        result.records, result.source_context, failed=result.status == "failed"
    )
    demand(
        result.counts == expected
        and all(type(v) is type(expected[k]) for k, v in result.counts.items()),
        "Inconsistent annotation counts.",
    )
    audit = {
        "outcome_reason": "incomplete"
        if result.status == "failed"
        else "no_identified_annotations"
        if not result.records
        else None,
        "targets": [],
        "activities": [],
        "assays": [],
    }
    demand(result.selection_audit == audit, "Invalid annotation audit.")
