"""Validate canonical UniProt accessions."""

import re

PATTERN = r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}"


def normalize_accession(value: str) -> str:
    accession = value.strip().upper()
    if not re.fullmatch(PATTERN, accession):
        raise ValueError(
            "Provide one canonical UniProt accession, such as P00699 (no isoform suffix)."
        )
    return accession
