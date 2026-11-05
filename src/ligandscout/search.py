"""Simple-mode source orchestration."""

from ligandscout.accessions import normalize_accession
from ligandscout.bindingdb import fetch_bindingdb


def run_simple(accession, *, client, sleep, now):
    return fetch_bindingdb(
        normalize_accession(accession), client=client, sleep=sleep, now=now
    )
