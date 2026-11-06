"""Reviewed UniProt ligand/cofactor annotations with sequence-site provenance."""

import hashlib
import logging
import re
from dataclasses import dataclass
from datetime import date

from ligandscout.accessions import normalize_accession
from ligandscout.api_http import Audit, request_json
from ligandscout.models import Measurement, SearchOutcome, SourceError, SourceResult
from ligandscout.uniprot_schema import UniProtSchemaError, annotation_counts

LOG = logging.getLogger(__name__)
GENERIC = {"substrate", "ligand", "unknown", "unspecified"}
CONTEXT = {
    "CATALYTIC ACTIVITY",
    "FUNCTION",
    "SUBUNIT",
    "ACTIVITY REGULATION",
    "BIOPHYSICOCHEMICAL PROPERTIES",
}


def obj(value):
    if not isinstance(value, dict):
        raise UniProtSchemaError("Expected an annotation object.")
    return value


def array(value):
    if not isinstance(value, list):
        raise UniProtSchemaError("Expected an annotation list.")
    return value


def text(value):
    if value is not None and not isinstance(value, str):
        raise UniProtSchemaError("Expected text or null.")
    return value if value and value.strip() else None


def chebi(value):
    value = text(value)
    if value is None:
        return None
    match = re.fullmatch(r"(?:ChEBI:)?CHEBI:([1-9][0-9]*)", value)
    if not match:
        raise UniProtSchemaError("Invalid ChEBI identifier.")
    return "CHEBI:" + match[1]


def evidence(values, path):
    result = []
    for index, item in enumerate(array(values)):
        item = obj(item)
        code = text(item.get("evidenceCode"))
        if not code or not re.fullmatch(r"ECO:[0-9]{7}", code):
            raise UniProtSchemaError("Invalid evidence code.")
        result.append(
            {
                "annotation_path": f"{path}/{index}",
                "evidence_code": code,
                "source": text(item.get("source")),
                "source_id": text(item.get("id")),
                "category": {
                    "ECO:0000269": "experimental",
                    "ECO:0000250": "similarity",
                    "ECO:0000255": "sequence_analysis",
                }.get(code, "other"),
            }
        )
    return result


def warn(warnings, code, path, message):
    warnings.append({"code": code, "message": f"{path}: {message}", "role": path})


def endpoint(value, warnings, path):
    value = obj(value)
    number = value.get("value")
    if number is not None and (type(number) is not int or number < 1):
        raise UniProtSchemaError("Invalid sequence coordinate.")
    if "modifier" in value and not text(value["modifier"]):
        raise UniProtSchemaError("Coordinate modifier must be nonempty text.")
    modifier = text(value.get("modifier")) or (
        "EXACT" if number is not None else "UNKNOWN"
    )
    if modifier == "EXACT" and number is None:
        raise UniProtSchemaError("Exact sequence coordinate has no position.")
    if modifier not in ("EXACT", "UNKNOWN", "UNSURE", "OUTSIDE", "BEFORE", "AFTER"):
        warn(
            warnings,
            "uncertain_coordinate",
            path,
            "Unrecognized coordinate modifier retained as uncertain.",
        )
    return {"value": number, "modifier": modifier}


def identity(ligand, warnings, path):
    ligand = obj(ligand)
    name, raw = text(ligand.get("name")), text(ligand.get("id"))
    chemical = None
    if raw:
        if raw.lower().startswith("chebi:"):
            chemical = chebi(raw)
        else:
            warn(
                warnings,
                "unsupported_ligand_identifier",
                path,
                "Unknown identifier namespace; no structure lookup.",
            )
    return chemical, name


def site(feature, path, primary, sequence, products, warnings):
    location = obj(feature.get("location", {}))
    start = endpoint(location.get("start", {}), warnings, path)
    end = endpoint(location.get("end", {}), warnings, path)
    seq_id = text(location.get("sequence"))
    canonical = seq_id in (None, primary)
    exact = start["modifier"] == end["modifier"] == "EXACT"
    if exact and (
        start["value"] > end["value"] or canonical and end["value"] > len(sequence)
    ):
        raise UniProtSchemaError(
            "Binding-site interval is outside its sequence or reversed."
        )
    if not exact:
        warn(
            warnings,
            "uncertain_coordinate",
            path,
            "Site coordinates are not exact; no residue sequence inferred.",
        )
    if not canonical:
        warn(
            warnings,
            "alternate_sequence",
            path,
            "Site belongs to an alternate sequence; canonical residues are not assigned.",
        )
    part = feature.get("ligandPart")
    if part is not None:
        part = obj(part)
        part = {"name": text(part.get("name")), "id": text(part.get("id"))}
        if part["id"] and part["id"].lower().startswith("chebi:"):
            part["id"] = chebi(part["id"])
    ligand = obj(feature.get("ligand", {}))
    for key in ("label", "note"):
        text(ligand.get(key))
    for crossref in array(feature.get("featureCrossReferences", [])):
        crossref = obj(crossref)
        if crossref.get("database") == "ChEBI":
            cid = chebi(crossref.get("id"))
            allowed = {ligand.get("id"), part["id"] if part else None}
            allowed = {
                chebi(x) for x in allowed if x and x.lower().startswith("chebi:")
            }
            if cid not in allowed:
                warn(
                    warnings,
                    "unassigned_feature_cross_reference",
                    path,
                    "Additional ChEBI reference retained without ligand assignment.",
                )
    return {
        "annotation_path": path,
        "feature_id": text(feature.get("featureId")),
        "location": {"start": start, "end": end, "sequence": seq_id},
        "sequence_scope": "canonical" if canonical else "alternate",
        "residue_sequence": sequence[start["value"] - 1 : end["value"]]
        if canonical and exact
        else None,
        "description": text(feature.get("description")),
        "ligand_part": part,
        "ligand_label": text(ligand.get("label")),
        "processed_products": [
            p
            for p in products
            if p["start"] <= start["value"] and p["end"] >= end["value"]
        ]
        if canonical and exact
        else [],
        "evidence": evidence(feature.get("evidences", []), path + "/evidences"),
    }


@dataclass
class ParsedEntry:
    records: list
    source_context: dict
    counts: dict
    selection_audit: dict
    warnings: list
    status: str


def parse_entry(data, *, query):
    data = obj(data)
    try:
        primary = normalize_accession(text(obj(data).get("primaryAccession")) or "")
        secondaries = [
            normalize_accession(x)
            for x in array(data.get("secondaryAccessions", []))
            if isinstance(x, str)
        ]
    except ValueError as exc:
        raise UniProtSchemaError("Invalid entry accession.") from exc
    if any(not isinstance(x, str) for x in data.get("secondaryAccessions", [])):
        raise UniProtSchemaError("Invalid secondary accession.")
    if primary != query and query not in secondaries:
        raise UniProtSchemaError("UniProt returned a different accession.")
    if data.get("entryType") != "UniProtKB reviewed (Swiss-Prot)":
        raise UniProtSchemaError("Only reviewed Swiss-Prot entries are supported.")
    seq = obj(data.get("sequence"))
    sequence = text(seq.get("value"))
    if (
        not sequence
        or not re.fullmatch(r"[A-Z]+", sequence)
        or type(seq.get("length")) is not int
        or seq["length"] != len(sequence)
    ):
        raise UniProtSchemaError("Invalid entry sequence or length.")
    organism = obj(data.get("organism"))
    scientific = text(organism.get("scientificName"))
    if (
        not scientific
        or type(organism.get("taxonId")) is not int
        or organism["taxonId"] < 1
    ):
        raise UniProtSchemaError("Invalid entry organism.")
    warnings = []
    try:
        name = text(
            obj(
                obj(
                    obj(data.get("proteinDescription", {})).get("recommendedName", {})
                ).get("fullName", {})
            ).get("value")
        )
    except UniProtSchemaError:
        name = None
    if not name:
        name = (
            data.get("uniProtkbId")
            if isinstance(data.get("uniProtkbId"), str)
            else primary
        )
        warn(
            warnings,
            "missing_protein_name",
            "/proteinDescription",
            "Using entry identifier as display name.",
        )
    audit = data.get("entryAudit", {})
    if not isinstance(audit, dict):
        audit = {}
        warn(
            warnings,
            "invalid_optional_metadata",
            "/entryAudit",
            "Entry audit unavailable.",
        )
    normalized_audit = {}
    for key in (
        "entryVersion",
        "sequenceVersion",
        "firstPublicDate",
        "lastAnnotationUpdateDate",
        "lastSequenceUpdateDate",
    ):
        value = audit.get(key)
        valid = value is None
        if key.endswith("Version"):
            valid = valid or type(value) is int and value > 0
        elif isinstance(value, str):
            try:
                date.fromisoformat(value)
                valid = True
            except ValueError:
                pass
        if not valid:
            warn(
                warnings,
                "invalid_optional_metadata",
                "/entryAudit/" + key,
                "Invalid optional audit value omitted.",
            )
        normalized_audit[key] = value if valid else None
    audit = normalized_audit
    if "references" in data and not isinstance(data["references"], list):
        warn(
            warnings,
            "invalid_optional_metadata",
            "/references",
            "Invalid reference list omitted.",
        )
    entry = {
        "primary_accession": primary,
        "query_accession": query,
        "secondary_accessions": secondaries,
        "entry_type": data["entryType"],
        "entry_name": data.get("uniProtkbId")
        if isinstance(data.get("uniProtkbId"), str)
        else None,
        "protein_name": name,
        "organism": scientific,
        "taxon_id": organism["taxonId"],
        "sequence": sequence,
        "sequence_length": len(sequence),
        "entry_audit": audit,
        "references": data.get("references", [])
        if isinstance(data.get("references", []), list)
        else [],
    }
    features, comments = (
        array(data.get("features", [])),
        array(data.get("comments", [])),
    )
    for f in features:
        if not isinstance(obj(f).get("type"), str):
            raise UniProtSchemaError("Feature type missing.")
    for c in comments:
        if not isinstance(obj(c).get("commentType"), str):
            raise UniProtSchemaError("Comment type missing.")
    products = []
    for i, feature in enumerate(features):
        if feature["type"] in ("Chain", "Peptide"):
            try:
                loc = obj(feature.get("location", {}))
                start = endpoint(loc.get("start", {}), warnings, f"/features/{i}")
                end = endpoint(loc.get("end", {}), warnings, f"/features/{i}")
                if (
                    loc.get("sequence") not in (None, primary)
                    or start["modifier"] != "EXACT"
                    or end["modifier"] != "EXACT"
                ):
                    continue
                if not 1 <= start["value"] <= end["value"] <= len(sequence):
                    raise UniProtSchemaError("Invalid processed product.")
                products.append(
                    {
                        "annotation_path": f"/features/{i}",
                        "name": text(feature.get("description")),
                        "start": start["value"],
                        "end": end["value"],
                    }
                )
            except UniProtSchemaError:
                warn(
                    warnings,
                    "invalid_optional_context",
                    f"/features/{i}",
                    "Processed product unavailable.",
                )
    context = {
        "entry": entry,
        "generic_sites": [],
        "generic_cofactors": [],
        "context_comments": [],
        "processed_products": products,
    }
    groups = {}

    def group(chemical, ligand_name, role):
        normalized = " ".join((ligand_name or "").split()).casefold()
        if not chemical and (not normalized or normalized in GENERIC):
            return None
        key = (
            chemical
            or f"UNIPROT:{primary}:NAME:{hashlib.sha256(normalized.encode()).hexdigest()}"
        )
        if key not in groups:
            groups[key] = Measurement(
                compound_id=key,
                compound_name=ligand_name,
                smiles=None,
                inchikey=None,
                target_description=name,
                evidence_type="uniprot_annotation",
                measurement_type=None,
                measurement_value=None,
                measurement_unit=None,
                pmid=None,
                doi=None,
                source_url=f"https://www.uniprot.org/uniprotkb/{primary}/entry",
                target_id=primary,
                target_accessions=[primary],
                target_organism=scientific,
                target_tax_id=str(organism["taxonId"]),
                target_match_status="entry_annotation",
                source_values={
                    "aliases": [],
                    "binding_annotations": [],
                    "cofactor_annotations": [],
                },
            )
        record = groups[key]
        if ligand_name and ligand_name not in record.source_values["aliases"]:
            record.source_values["aliases"].append(ligand_name)
            record.compound_name = record.compound_name or ligand_name
        if role not in record.ligand_roles:
            record.ligand_roles.append(role)
        return record

    for i, feature in enumerate(features):
        if feature["type"] != "Binding site":
            continue
        path = f"/features/{i}"
        chemical, ligand_name = identity(feature.get("ligand", {}), warnings, path)
        binding = site(feature, path, primary, sequence, products, warnings)
        record = group(chemical, ligand_name, "binding_site")
        if record is None:
            context["generic_sites"].append(binding)
        else:
            record.binding_sites.append(binding)
            record.annotation_evidence.extend(binding["evidence"])
            record.source_values["binding_annotations"].append(
                {"annotation_path": path, "data": feature}
            )
    for i, comment in enumerate(comments):
        path = f"/comments/{i}"
        kind = comment["commentType"]
        if kind in CONTEXT or kind == "COFACTOR":
            try:
                if "texts" in comment:
                    for j, t in enumerate(array(comment["texts"])):
                        text(obj(t).get("value"))
                        evidence(t.get("evidences", []), f"{path}/texts/{j}/evidences")
                if kind == "CATALYTIC ACTIVITY":
                    reaction = obj(comment.get("reaction"))
                    text(reaction.get("name"))
                    text(reaction.get("ecNumber"))
                    evidence(
                        reaction.get("evidences", []), path + "/reaction/evidences"
                    )
                    for crossref in array(reaction.get("reactionCrossReferences", [])):
                        crossref = obj(crossref)
                        text(crossref.get("database"))
                        text(crossref.get("id"))
                context["context_comments"].append(
                    {"comment_type": kind, "annotation_path": path, "data": comment}
                )
            except UniProtSchemaError:
                warn(
                    warnings,
                    "invalid_optional_context",
                    path,
                    "Context comment unavailable.",
                )
        if kind != "COFACTOR":
            continue
        if "note" in comment:
            for note in array(obj(comment["note"]).get("texts", [])):
                note = obj(note)
                text(note.get("value"))
                evidence(note.get("evidences", []), path + "/note/texts/evidences")
        parent_evidence = evidence(comment.get("evidences", []), path + "/evidences")
        for j, item in enumerate(array(comment.get("cofactors", []))):
            item = obj(item)
            item_path = f"{path}/cofactors/{j}"
            cross = obj(item.get("cofactorCrossReference", {}))
            chemical = (
                chebi(cross.get("id")) if cross.get("database") == "ChEBI" else None
            )
            ligand_name = text(item.get("name"))
            item_evidence = (
                evidence(item.get("evidences", []), item_path + "/evidences")
                + parent_evidence
            )
            original = {
                "annotation_path": item_path,
                "data": item,
                "note": comment.get("note"),
                "evidence": item_evidence,
            }
            record = group(chemical, ligand_name, "cofactor")
            if record is None:
                context["generic_cofactors"].append(original)
            else:
                record.annotation_evidence.extend(item_evidence)
                record.source_values["cofactor_annotations"].append(original)
    records = sorted(groups.values(), key=lambda r: r.compound_id)
    for record in records:
        record.ligand_roles.sort(key=("binding_site", "cofactor").index)
        if any(
            item["category"] in ("similarity", "sequence_analysis")
            for item in record.annotation_evidence
        ):
            warn(
                warnings,
                "inferred_annotation",
                record.compound_id,
                "Includes similarity/sequence-based evidence; not all annotations are experimental.",
            )
    status = "success" if records else "no_matches"
    return ParsedEntry(
        records,
        context,
        annotation_counts(records, context),
        {
            "outcome_reason": None if records else "no_identified_annotations",
            "targets": [],
            "activities": [],
            "assays": [],
        },
        warnings,
        status,
    )


def fetch_uniprot(query, *, client, sleep, now, monotonic):
    from ligandscout.ccd import enrich_ligands
    from ligandscout.go_annotations import enrich_go

    LOG.info("Searching UniProt for %s", query)
    audit = Audit()
    url = f"https://rest.uniprot.org/uniprotkb/{query}.json"
    response = request_json(
        url,
        params={},
        role="entry",
        chemical_id=None,
        audit=audit,
        client=client,
        sleep=sleep,
        now=now,
        monotonic=monotonic,
    )
    result = SourceResult(
        schema_version=4,
        accession=query,
        source="UniProt",
        retrieved_at=response.retrieved_at,
        request={"url": url, "parameters": {}},
        source_version=response.headers.get("x-uniprot-release"),
        status="failed",
        records=[],
        error=response.error,
        policy="uniprot-annotated-ligands-v1",
        counts=annotation_counts([], {}, failed=True),
        selection_audit={
            "outcome_reason": "incomplete",
            "targets": [],
            "activities": [],
            "assays": [],
        },
    )
    if response.error is None:
        try:
            parsed = parse_entry(response.data, query=query)
        except (UniProtSchemaError, ValueError) as exc:
            result.error = SourceError(
                "schema", f"UniProt entry validation failed: {exc}", response.status
            )
        else:
            result.records, result.source_context = (
                parsed.records,
                parsed.source_context,
            )
            result.counts, result.selection_audit = (
                parsed.counts,
                parsed.selection_audit,
            )
            result.warnings, result.status = parsed.warnings, parsed.status
            enrich_go(
                result,
                response.data,
                audit=audit,
                client=client,
                sleep=sleep,
                now=now,
                monotonic=monotonic,
            )
            enrich_ligands(
                result,
                audit=audit,
                client=client,
                sleep=sleep,
                now=now,
                monotonic=monotonic,
            )
    LOG.info(
        "UniProt %s: %s; %s ligand annotation records",
        query,
        result.status,
        len(result.records),
    )
    return SearchOutcome(result, response.body, audit.exchanges)
