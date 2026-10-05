"""Shared fixtures. No test may reach the network, and none needs real Portal credentials."""

import socket

import pytest

from igvf_portal.connection import PConnection
from igvf_portal.enums import IgvfMode
from igvf_portal.types import IgvfRecord


@pytest.fixture(autouse=True)
def no_network(monkeypatch: pytest.MonkeyPatch) -> None:
    """Fail any test that tries to open a network connection (e.g. to the IGVF Portal)."""

    def _refuse(*_args: object, **_kwargs: object) -> None:
        raise RuntimeError("Network access is disabled in tests.")

    monkeypatch.setattr(socket.socket, "connect", _refuse)
    monkeypatch.setattr(socket, "create_connection", _refuse)


@pytest.fixture(autouse=True)
def fake_credentials(monkeypatch: pytest.MonkeyPatch) -> None:
    """Set dummy Portal credentials, so code that checks for them can run."""
    monkeypatch.setenv("IGVF_API_KEY", "test-api-key")
    monkeypatch.setenv("IGVF_SECRET_KEY", "test-secret-key")


@pytest.fixture
def connection() -> PConnection:
    """A PConnection to prod. Constructing it makes no requests; seed records with seed_record."""
    return PConnection.new(IgvfMode.prod)


def seed_record(connection: PConnection, record: IgvfRecord, frame: str | None = None) -> None:
    """Put a record in the connection's lookup cache, so lookup_record finds it offline.

    The record is cached under its accession, @id, and aliases.
    """
    connection._cache_record(
        keys=[record["@id"]],
        frame=frame,
        database=connection.submission,
        record=record,
    )


class FakeProperty:
    """Just enough of an igvf_utils schema property: its name and JSON schema."""

    def __init__(self, name: str, schema: dict[str, object]) -> None:
        """Make a property with the given name and JSON schema."""
        self.name = name
        self.schema = schema


class FakeSchema:
    """Just enough of igvf_utils' IgvfSchema for building register payloads."""

    name = "tabular_file"

    def __init__(self, properties: dict[str, dict[str, object]]) -> None:
        """Make a schema from a map of property name to JSON schema."""
        self._properties = {name: FakeProperty(name, s) for name, s in properties.items()}

    @property
    def properties(self) -> list[FakeProperty]:
        """All properties of the schema."""
        return list(self._properties.values())

    def get_property_from_name(self, name: str) -> FakeProperty:
        """Get the named property."""
        return self._properties[name]


FAKE_SCHEMA = FakeSchema(
    {
        "aliases": {"type": "array", "items": {"type": "string"}},
        "file_size": {"type": "integer"},
        "controlled_access": {"type": "boolean"},
        "description": {"type": "string"},
        "cell_qualifier": {"type": "string"},
        "file_set": {"type": "string"},
        "attachment": {"type": "object"},
        "derived_from": {"type": "array", "items": {"type": "string"}},
    }
)
