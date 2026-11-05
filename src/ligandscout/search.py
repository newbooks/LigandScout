"""Simple-mode source orchestration."""

from dataclasses import replace
from time import monotonic as system_monotonic

from ligandscout.accessions import normalize_accession
from ligandscout.bindingdb import fetch_bindingdb
from ligandscout.chembl import fetch_chembl
from ligandscout.models import RunOutcome, upgrade_bindingdb


def run_simple(accession, *, client, sleep, now, monotonic=system_monotonic):
    accession = normalize_accession(accession)
    bindingdb = fetch_bindingdb(accession, client=client, sleep=sleep, now=now)
    bindingdb = replace(bindingdb, result=upgrade_bindingdb(bindingdb.result))
    chembl = fetch_chembl(
        accession, client=client, sleep=sleep, now=now, monotonic=monotonic
    )
    return RunOutcome(accession, [bindingdb, chembl])
