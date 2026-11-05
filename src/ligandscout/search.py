"""Search orchestration for BindingDB and UniProt."""

from dataclasses import replace
from time import monotonic as system_monotonic

from ligandscout.accessions import normalize_accession
from ligandscout.bindingdb import fetch_bindingdb
from ligandscout.models import RunOutcome, upgrade_result
from ligandscout.uniprot import fetch_uniprot


def run_search(accession, *, client, sleep, now, monotonic=system_monotonic):
    accession = normalize_accession(accession)
    outcome = fetch_bindingdb(accession, client=client, sleep=sleep, now=now)
    uniprot = fetch_uniprot(
        accession, client=client, sleep=sleep, now=now, monotonic=monotonic
    )
    return RunOutcome(
        accession, [replace(outcome, result=upgrade_result(outcome.result)), uniprot]
    )
