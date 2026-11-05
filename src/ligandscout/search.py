"""Simple-mode orchestration for the implemented BindingDB adapter."""

from dataclasses import replace

from ligandscout.accessions import normalize_accession
from ligandscout.bindingdb import fetch_bindingdb
from ligandscout.models import RunOutcome, upgrade_result


def run_simple(accession, *, client, sleep, now):
    accession = normalize_accession(accession)
    outcome = fetch_bindingdb(accession, client=client, sleep=sleep, now=now)
    return RunOutcome(
        accession, [replace(outcome, result=upgrade_result(outcome.result))]
    )
