"""Exact ChEBI-to-CCD structure enrichment; no interaction inference."""

import json
import re
from urllib.parse import quote

from ligandscout.api_http import request_json
from ligandscout.uniprot_schema import annotation_counts

SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
COMPONENT = "https://data.rcsb.org/rest/v1/core/chemcomp/"


class EnrichmentError(ValueError):
    """Incomplete or malformed chemical enrichment."""


def lookup(chemical, *, audit, client, sleep, now, monotonic):
    request_ids = []

    def request(url, params, role):
        response = request_json(
            url,
            params=params,
            role=role,
            chemical_id=chemical,
            audit=audit,
            client=client,
            sleep=sleep,
            now=now,
            monotonic=monotonic,
        )
        request_ids.append(audit.request_count)
        if response.error:
            raise EnrichmentError(response.error.message)
        return response

    provenance = {
        "provider": "RCSB CCD",
        "status": "unavailable",
        "requested_chemical_id": chemical,
        "smiles_field": None,
        "candidate_components": [],
        "request_ids": request_ids,
    }
    try:
        identifiers, total, start = [], None, 0
        while True:
            query = {
                "query": {
                    "type": "terminal",
                    "service": "text_chem",
                    "parameters": {
                        "attribute": "rcsb_chem_comp_related.resource_accession_code",
                        "operator": "exact_match",
                        "value": chemical,
                    },
                },
                "return_type": "mol_definition",
                "request_options": {"paginate": {"start": start, "rows": 100}},
            }
            response = request(
                SEARCH, {"json": json.dumps(query, separators=(",", ":"))}, "ccd_search"
            )
            if response.status == 204:
                if start:
                    raise EnrichmentError("CCD search became empty during pagination.")
                provenance["status"] = "no_exact_mapping"
                return provenance, None, None, None
            data = response.data
            if (
                not isinstance(data, dict)
                or data.get("result_type") != "mol_definition"
                or type(data.get("total_count")) is not int
                or data["total_count"] < 0
            ):
                raise EnrichmentError("Invalid CCD search response.")
            if total is not None and total != data["total_count"]:
                raise EnrichmentError("CCD search total changed.")
            total = data["total_count"]
            rows = data.get("result_set", [])
            if not isinstance(rows, list):
                raise EnrichmentError("Invalid CCD result set.")
            page = []
            for row in rows:
                identifier = row.get("identifier") if isinstance(row, dict) else None
                if not isinstance(identifier, str) or not re.fullmatch(
                    r"[A-Za-z0-9]+", identifier
                ):
                    raise EnrichmentError("Invalid CCD component ID.")
                if identifier in identifiers or identifier in page:
                    raise EnrichmentError("Duplicate CCD component during pagination.")
                page.append(identifier)
            identifiers.extend(page)
            if (
                len(identifiers) > total
                or len(page) > 100
                or (len(identifiers) < total and len(page) != 100)
            ):
                raise EnrichmentError("Incomplete CCD search page.")
            if len(identifiers) == total:
                break
            start += len(page)
        eligible = []
        for identifier in identifiers:
            response = request(
                COMPONENT + quote(identifier, safe=""), {}, "ccd_component"
            )
            data = response.data
            if (
                not isinstance(data, dict)
                or data.get("rcsb_id") != identifier
                or not isinstance(data.get("chem_comp"), dict)
                or data["chem_comp"].get("id") != identifier
            ):
                raise EnrichmentError("CCD component identity mismatch.")
            related = data.get("rcsb_chem_comp_related", [])
            if not isinstance(related, list) or any(
                not isinstance(r, dict)
                or not isinstance(r.get("resource_name"), str)
                or not isinstance(r.get("resource_accession_code"), str)
                for r in related
            ):
                raise EnrichmentError("Invalid CCD chemical cross-references.")
            matches = [
                r
                for r in related
                if r["resource_name"] == "ChEBI"
                and r["resource_accession_code"] == chemical
            ]
            if not matches:
                provenance["candidate_components"].append(
                    {
                        "ccd_id": identifier,
                        "mapping": None,
                        "smiles_field": None,
                        "smiles": None,
                        "inchikey": None,
                        "request_id": audit.request_count,
                        "reason": "no_exact_mapping",
                    }
                )
                continue
            field, smiles, inchikey = descriptors(data)
            candidate = {
                "ccd_id": identifier,
                "mapping": matches[0],
                "smiles_field": field,
                "smiles": smiles,
                "inchikey": inchikey,
                "request_id": audit.request_count,
                "reason": None,
            }
            eligible.append(candidate)
            provenance["candidate_components"].append(candidate)
        if not eligible:
            provenance["status"] = "no_exact_mapping"
        elif any(c["smiles"] is None for c in eligible):
            provenance["status"] = "no_descriptor"
        elif (
            len({(c["smiles_field"], c["smiles"], c["inchikey"]) for c in eligible})
            != 1
        ):
            provenance["status"] = "ambiguous"
        else:
            chosen = eligible[0]
            provenance.update(status="available", smiles_field=chosen["smiles_field"])
            return provenance, chosen["smiles"], chosen["inchikey"], None
        return provenance, None, None, None
    except EnrichmentError as exc:
        provenance["status"] = "unavailable"
        return provenance, None, None, str(exc)


def enrich_ligands(result, *, audit, client, sleep, now, monotonic):
    cache = {}
    for record in result.records:
        chemical = (
            record.compound_id if record.compound_id.startswith("CHEBI:") else None
        )
        if chemical:
            if chemical not in cache:
                cache[chemical] = lookup(
                    chemical,
                    audit=audit,
                    client=client,
                    sleep=sleep,
                    now=now,
                    monotonic=monotonic,
                )
            provenance, record.smiles, record.inchikey, error = cache[chemical]
        else:
            provenance = {
                "provider": "RCSB CCD",
                "status": "not_applicable",
                "requested_chemical_id": None,
                "smiles_field": None,
                "candidate_components": [],
                "request_ids": [],
            }
            error = None
        record.structure_provenance = provenance
        if provenance["status"] not in ("available", "not_applicable"):
            result.warnings.append(
                {
                    "code": "structure_" + provenance["status"],
                    "message": f"{record.compound_id}: SMILES {provenance['status']}."
                    + (f" {error}" if error else ""),
                    "role": record.compound_id,
                }
            )
    result.counts = annotation_counts(result.records, result.source_context)
    return result


def descriptors(data):
    """Select a descriptor from an already identified CCD component."""
    values = data.get("rcsb_chem_comp_descriptor", {})
    if not isinstance(values, dict) or any(
        v is not None and not isinstance(v, str)
        for v in (
            values.get("SMILES_stereo"),
            values.get("SMILES"),
            values.get("InChIKey"),
        )
    ):
        raise EnrichmentError("Invalid CCD descriptor.")
    field = next(
        (k for k in ("SMILES_stereo", "SMILES") if (values.get(k) or "").strip()), None
    )
    return (
        field,
        values[field] if field else None,
        values.get("InChIKey") if (values.get("InChIKey") or "").strip() else None,
    )
