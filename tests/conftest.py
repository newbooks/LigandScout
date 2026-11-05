"""Deterministic BindingDB fixtures."""

import json
from datetime import datetime, timezone
from pathlib import Path

import httpx
import pytest

from ligandscout.bindingdb import fetch_bindingdb


@pytest.fixture
def now():
    return lambda: datetime(2026, 11, 4, 15, tzinfo=timezone.utc)


@pytest.fixture
def payload():
    return json.loads(
        (Path(__file__).parent / "fixtures/bindingdb-positive.json").read_text()
    )


@pytest.fixture
def fetch(now):
    def call(handler):
        sleeps = []
        with httpx.Client(transport=httpx.MockTransport(handler)) as client:
            outcome = fetch_bindingdb(
                "P00699", client=client, sleep=sleeps.append, now=now
            )
            assert not client.is_closed
        return outcome, sleeps

    return call


@pytest.fixture
def positive(fetch, payload):
    return fetch(lambda request: httpx.Response(200, json=payload))[0]


@pytest.fixture
def uniprot_entry():
    return json.loads(
        (Path(__file__).parent / "fixtures/uniprot-annotations.json").read_text()
    )


@pytest.fixture
def empty_uniprot(uniprot_entry):
    uniprot_entry.update(primaryAccession="P00699", features=[], comments=[])
    return uniprot_entry
