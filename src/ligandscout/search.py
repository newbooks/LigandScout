"""Search orchestration for BindingDB, UniProt and RCSB PDB."""

from dataclasses import replace
from time import monotonic as system_monotonic

from ligandscout.accessions import normalize_accession
from ligandscout.bindingdb import fetch_bindingdb
from ligandscout.cache import load_completed
from ligandscout.models import RunOutcome, upgrade_result
from ligandscout.rcsb import fetch_rcsb
from ligandscout.uniprot import fetch_uniprot


def run_search(
    accession, *, client, sleep, now, monotonic=system_monotonic, root=None, force=False
):
    accession = normalize_accession(accession)

    def cached(source):
        return (
            load_completed(root, accession, source)
            if root is not None and not force
            else None
        )

    outcome = cached("BindingDB")
    if outcome is None:
        outcome = fetch_bindingdb(accession, client=client, sleep=sleep, now=now)
    uniprot = cached("UniProt")
    if uniprot is None:
        uniprot = fetch_uniprot(
            accession, client=client, sleep=sleep, now=now, monotonic=monotonic
        )
    rcsb = cached("RCSB PDB")
    if rcsb is None:
        rcsb = fetch_rcsb(
            accession,
            client=client,
            sleep=sleep,
            now=now,
            monotonic=monotonic,
            uniprot_context=uniprot.result.source_context.get("entry")
            if uniprot.result.error is None
            else None,
        )
    return RunOutcome(
        accession,
        [replace(s, result=upgrade_result(s.result)) for s in (outcome, uniprot, rcsb)],
    )
