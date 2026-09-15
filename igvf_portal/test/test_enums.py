from typing import cast

import pytest

from igvf_portal.enums import ContentType
from igvf_portal.types import IgvfRecord


def _record(**fields: str) -> IgvfRecord:
    return cast(
        IgvfRecord,
        {
            "status": "released",
            "upload_status": "validated",
            "content_type": "fragments",
            "file_format": "bed",
        }
        | fields,
    )


@pytest.mark.parametrize(
    ("fields", "wanted"),
    [
        ({}, True),
        ({"file_format": "tsv"}, True),
        ({"file_format": "csv"}, False),
        ({"content_type": "peaks"}, False),
        ({"status": "deleted"}, False),
        ({"upload_status": "invalidated"}, False),
        ({"upload_status": "validation exempted"}, True),
    ],
)
def test_is_wanted(fields: dict[str, str], wanted: bool) -> None:
    assert ContentType.FRAGMENTS.is_wanted(_record(**fields)) is wanted


def test_is_wanted_without_upload_status() -> None:
    record = _record()
    del record["upload_status"]
    assert ContentType.FRAGMENTS.is_wanted(record)


def test_every_content_type_has_formats_and_extension() -> None:
    for content_type in ContentType:
        assert len(content_type.file_formats) > 0
        assert content_type.extension
