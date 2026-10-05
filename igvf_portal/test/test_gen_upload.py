"""Tests for generating the upload TSVs and script (gen-upload-script, used by IGVF_UPLOAD)."""

import json
import logging
from pathlib import Path
from typing import cast

import pytest

from igvf_portal import utils
from igvf_portal.connection import PConnection
from igvf_portal.enums import OutputCategory
from igvf_portal.gen_upload_config import GenUploadConfig
from igvf_portal.igvf_uploads import IgvfDocument, IgvfPseudobulk
from igvf_portal.parallel_logger import ParallelLogger
from igvf_portal.tools import register
from igvf_portal.types import IgvfRecord, PseudobulkId, UploadRow
from igvf_portal.upload_state import UploadState

from .conftest import FAKE_SCHEMA, seed_record

_ANNOTATION_COLUMNS = (
    "pseudobulk_id",
    "cell_name",
    "cleaned_cell_name",
    "CL_id",
    "cell_description",
    "CL_term_name",
    "subsample",
    "cleaned_subsample",
)


def _write_annotations(path: Path, rows: list[dict[str, str]]) -> Path:
    columns = list(rows[0])
    lines = ["\t".join(columns)] + ["\t".join(row[c] for c in columns) for row in rows]
    path.write_text("\n".join(lines) + "\n")
    return path


def _annotation(pseudobulk_id: str, **fields: str) -> dict[str, str]:
    row = {
        "pseudobulk_id": pseudobulk_id,
        "cell_name": "naive T cell",
        "cleaned_cell_name": "naive_T_cell",
        "CL_id": "CL:0000084",
        "cell_description": "T cell",
        "CL_term_name": "T cell",
        "subsample": "IGVFSM0000SSSS",
        "cleaned_subsample": "IGVFSM0000SSSS",
    }
    return row | fields


def _config(basedir: Path, connection: PConnection, annotations: Path) -> GenUploadConfig:
    return GenUploadConfig(
        basedir=basedir,
        input_file_sets="IGVFDS0000AAAA",
        connection=connection,
        annotations_path=annotations,
        logger=ParallelLogger.new(utils.get_logger("test-gen-upload", level=logging.INFO)),
    )


@pytest.fixture
def seeded_connection(connection: PConnection) -> PConnection:
    """A connection that knows the input analysis set and the T cell sample term."""
    seed_record(
        connection,
        cast(
            IgvfRecord,
            {
                "@id": "/analysis-sets/IGVFDS0000AAAA/",
                "accession": "IGVFDS0000AAAA",
                "aliases": ["other-lab:analysis"],
                "files": [],
            },
        ),
    )
    seed_record(
        connection,
        cast(
            IgvfRecord,
            {"@id": "/sample-terms/CL_0000084/", "aliases": [], "term_name": "T cell"},
        ),
        frame="object",
    )
    return connection


@pytest.mark.parametrize(
    "extra_columns",
    [
        {},
        # PSEUDOBULK_RNA writes cell_qualifier when the metadata has it
        {"cell_qualifier": "naive"},
        {"cell_qualifier": ""},
        {"unexpected": "x"},
    ],
)
def test_annotations(
    tmp_path: Path, connection: PConnection, extra_columns: dict[str, str]
) -> None:
    annotations = _write_annotations(
        tmp_path / "cell_name_to_annotation_mapping.tsv",
        [
            _annotation("naive_T_cell-a", **extra_columns),
            _annotation("B_cell-a", **extra_columns),
        ],
    )
    config = _config(tmp_path, connection, annotations)
    assert set(config.annotations) == {"naive_T_cell-a", "B_cell-a"}
    row = config.annotations[PseudobulkId("naive_T_cell-a")]
    expected_columns = set(_ANNOTATION_COLUMNS) | (
        {"cell_qualifier"} if "cell_qualifier" in extra_columns else set()
    )
    assert set(row) == expected_columns
    assert row.get("cell_qualifier") == extra_columns.get("cell_qualifier")


def test_annotations_missing_column(tmp_path: Path, connection: PConnection) -> None:
    row = _annotation("naive_T_cell-a")
    del row["CL_id"]
    annotations = _write_annotations(tmp_path / "annotations.tsv", [row])
    with pytest.raises(ValueError, match="missing required columns: CL_id"):
        _ = _config(tmp_path, connection, annotations).annotations


@pytest.mark.parametrize(
    ("cell_qualifier", "expected"),
    [
        (None, "naive"),  # inferred from cell_name minus the term name
        ("", "naive"),  # an empty cell is inferred too
        ("activated", "activated"),
    ],
)
def test_pseudobulk_row(
    tmp_path: Path,
    seeded_connection: PConnection,
    cell_qualifier: str | None,
    expected: str,
) -> None:
    extra = {} if cell_qualifier is None else {"cell_qualifier": cell_qualifier}
    annotations = _write_annotations(
        tmp_path / "annotations.tsv", [_annotation("naive_T_cell-a", **extra)]
    )
    folder = tmp_path / "pseudobulks" / "naive_T_cell-a"
    folder.mkdir(parents=True)
    config = _config(tmp_path, seeded_connection, annotations)
    row = IgvfPseudobulk().get_row(check_path=folder, config=config, doc_aliases=[])
    assert row == {
        "aliases": "anshul-kundaje:pseudobulk-IGVFDS0000AAAA-naive_T_cell-IGVFSM0000SSSS",
        "award": "/awards/HG012069/",
        "lab": "/labs/anshul-kundaje/",
        "file_set_type": "pseudobulk analysis",
        "cell_type": "/sample-terms/CL_0000084/",
        "samples": "IGVFSM0000SSSS",
        "input_file_sets": "/analysis-sets/IGVFDS0000AAAA/",
        "documents": "",
        "merged": False,
        "cell_qualifier": expected,
    }


def test_document_row_attachment(tmp_path: Path, seeded_connection: PConnection) -> None:
    """The attachment is JSON naming the document's path relative to the upload folder."""
    annotations = _write_annotations(tmp_path / "annotations.tsv", [_annotation("pb-a")])
    folder = tmp_path / "pseudobulks" / "pb-a"
    folder.mkdir(parents=True)
    (folder / "qc_plot.pdf").write_bytes(b"%PDF")
    document = IgvfDocument(
        output_category=OutputCategory.PSEUDOBULK,
        match_glob="*.pdf",
        document_type="plate map",
        description="QC plot",
    )
    doc_aliases = []
    row = document.get_row(
        check_path=folder,
        config=_config(tmp_path, seeded_connection, annotations),
        doc_aliases=doc_aliases,
    )
    assert row is not None
    assert json.loads(row["attachment"]) == {"path": "pseudobulks/pb-a/qc_plot.pdf"}
    assert doc_aliases == [row["aliases"]]


def test_written_tsv_round_trips_through_register(
    tmp_path: Path, seeded_connection: PConnection
) -> None:
    """The upload TSVs are read back by register: values with quotes must survive intact."""
    annotations = _write_annotations(tmp_path / "annotations.tsv", [_annotation("pb-a")])
    (tmp_path / "pseudobulks" / "pb-a").mkdir(parents=True)
    upload_state = UploadState(
        basedir=tmp_path, config=_config(tmp_path, seeded_connection, annotations)
    )
    rows = [
        cast(
            UploadRow,
            {
                "aliases": "lab:a",
                "award": "/awards/A/",
                "lab": "/labs/a/",
                "description": 'CD4 "naive" T cells, from 2 donors',
                "cell_qualifier": "it's naive",
            },
        )
    ]
    upload_state._write_tsv("tabular_file", rows, analysis_step=None)
    tsv = tmp_path / "upload_tsvs" / "tabular_file.0.tsv"
    payloads = list(
        register._iter_payloads(
            schema=cast(register.IgvfSchema, FAKE_SCHEMA),
            infile=tsv,
            drop_extra_fields=True,
        )
    )
    assert payloads == [
        {
            "_profile": "tabular_file",
            "aliases": ["lab:a"],
            "description": 'CD4 "naive" T cells, from 2 donors',
            "cell_qualifier": "it's naive",
        }
    ]
    assert upload_state.submission_rows[-1] == (
        'igvf-portal register $dry_run_arg --igvf-mode "$igvf_mode" --profile-id tabular_file'
        ' --infile "upload_tsvs/tabular_file.0.tsv" --log-level "info"'
    )
