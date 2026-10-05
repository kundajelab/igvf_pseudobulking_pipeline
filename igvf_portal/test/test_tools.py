"""Tests for the tools used by the Nextflow modules, run offline against seeded records."""

import io
import json
import logging
from pathlib import Path
from typing import Final, cast

import pytest

from igvf_portal.connection import PConnection
from igvf_portal.enums import ContentType, IgvfMode, MultipleRecordsAction
from igvf_portal.tools import find_bad_beds, get_species, get_url, register
from igvf_portal.types import AccessionId, IgvfRecord, PortalId

from .conftest import FAKE_SCHEMA, seed_record

_LOGGER = logging.getLogger("test")


_FILE_FORMATS: dict[str, str] = {
    content_type.value: next(iter(sorted(content_type.file_formats)))
    for content_type in ContentType
}
"""The format get-url downloads for each content type; other content types are given tsv."""


def _file(accession: str, content_type: str, **fields: object) -> IgvfRecord:
    """A file record as embedded in its file set, with the format get-url wants by default."""
    return cast(
        IgvfRecord,
        {
            "@id": f"/files/{accession}/",
            "@type": ["File", "Item"],
            "accession": accession,
            "aliases": [],
            "status": "released",
            "upload_status": "validated",
            "content_type": content_type,
            "file_format": _FILE_FORMATS.get(content_type, "tsv"),
            "file_set": {"accession": "IGVFDS0000AAAA"},
            "href": f"/files/{accession}/@@download/{accession}.h5ad",
            "s3_uri": f"s3://igvf-public/2025/01/01/uuid/{accession}.h5ad",
        }
        | fields,
    )


_PSEUDOBULK_SET: Final[dict[str, str]] = {"accession": "IGVFDS0000PPPP"}


def _analysis_set(*files: IgvfRecord, accession: str = "IGVFDS0000AAAA") -> IgvfRecord:
    return cast(
        IgvfRecord,
        {
            "@id": f"/analysis-sets/{accession}/",
            "@type": ["AnalysisSet", "FileSet", "Item"],
            "accession": accession,
            "aliases": [],
            "status": "released",
            "files": list(files),
        },
    )


def _content_type_record(
    connection: PConnection,
    analysis_set: IgvfRecord,
    action: MultipleRecordsAction = MultipleRecordsAction.KEEP_UNFILTERED,
) -> IgvfRecord | None:
    return get_url._get_content_type_record(
        analysis_set_record=analysis_set,
        content_type=ContentType.MATRIX,
        multiple_records_action=action,
        connection=connection,
        logger=_LOGGER,
    )


def test_get_content_type_record_single(connection: PConnection) -> None:
    matrix = _file("IGVFFI0000MMMM", "cell by gene matrix")
    fragments = _file("IGVFFI0000FFFF", "fragments")
    deleted = _file("IGVFFI0000DDDD", "cell by gene matrix", status="deleted")
    analysis_set = _analysis_set(matrix, fragments, deleted)
    assert _content_type_record(connection, analysis_set) is matrix


def test_get_content_type_record_skips_invalidated(connection: PConnection) -> None:
    """An invalidated file is not a candidate, but is no reason to fail either."""
    matrix = _file("IGVFFI0000MMMM", "cell by gene matrix")
    invalidated = _file("IGVFFI0000IIII", "cell by gene matrix", upload_status="invalidated")
    assert _content_type_record(connection, _analysis_set(invalidated, matrix)) is matrix
    assert _content_type_record(connection, _analysis_set(invalidated)) is None


@pytest.mark.parametrize("file_format", ["tar", "hdf5"])
def test_get_content_type_record_ignores_wrong_format(
    connection: PConnection, file_format: str
) -> None:
    """Analysis sets can hold other files: only matrices the pipeline can read are candidates."""
    other = _file("IGVFFI0000OOOO", "cell by gene matrix", file_format=file_format)
    matrix = _file("IGVFFI0000MMMM", "cell by gene matrix")
    assert _content_type_record(connection, _analysis_set(other, matrix)) is matrix
    assert _content_type_record(connection, _analysis_set(other)) is None


@pytest.mark.parametrize("file_format", ["tsv", "bed"])
def test_get_content_type_record_fragments_formats(
    connection: PConnection, file_format: str
) -> None:
    fragments = _file("IGVFFI0000FFFF", "fragments", file_format=file_format)
    record = get_url._get_content_type_record(
        analysis_set_record=_analysis_set(fragments),
        content_type=ContentType.FRAGMENTS,
        multiple_records_action=MultipleRecordsAction.KEEP_UNFILTERED,
        connection=connection,
        logger=_LOGGER,
    )
    assert record is fragments


def test_get_content_type_record_without_upload_status(connection: PConnection) -> None:
    """A file record without upload_status is still a candidate."""
    matrix = _file("IGVFFI0000MMMM", "cell by gene matrix")
    del matrix["upload_status"]
    assert _content_type_record(connection, _analysis_set(matrix)) is matrix


def test_get_content_type_record_none(connection: PConnection) -> None:
    analysis_set = _analysis_set(_file("IGVFFI0000FFFF", "fragments"))
    assert _content_type_record(connection, analysis_set) is None


@pytest.mark.parametrize(
    ("action", "expected"),
    [
        (MultipleRecordsAction.KEEP_UNFILTERED, "IGVFFI0000UUUU"),
        (MultipleRecordsAction.KEEP_FILTERED, "IGVFFI0000FFFF"),
    ],
)
def test_get_content_type_record_multiple(
    connection: PConnection, action: MultipleRecordsAction, expected: str
) -> None:
    """With several matrices, the filtered property (only in each file's own record) decides."""
    unfiltered = _file("IGVFFI0000UUUU", "cell by gene matrix")
    filtered = _file("IGVFFI0000FFFF", "cell by gene matrix")
    # the object frame has filtered, which the embedded file records lack
    for accession, is_filtered in (("IGVFFI0000UUUU", False), ("IGVFFI0000FFFF", True)):
        seed_record(
            connection,
            _file(accession, "cell by gene matrix", filtered=is_filtered),
            frame="object",
        )
    record = _content_type_record(connection, _analysis_set(unfiltered, filtered), action)
    assert record is not None
    assert record["accession"] == expected


def test_get_content_type_record_multiple_raise(connection: PConnection) -> None:
    analysis_set = _analysis_set(
        _file("IGVFFI0000UUUU", "cell by gene matrix"),
        _file("IGVFFI0000FFFF", "cell by gene matrix"),
    )
    with pytest.raises(ValueError, match="Found 2 cell by gene matrix files"):
        _content_type_record(connection, analysis_set, MultipleRecordsAction.RAISE)


def test_get_content_type_record_multiple_none_kept(connection: PConnection) -> None:
    """If every candidate is filtered out, raise rather than quietly download no matrix."""
    first = _file("IGVFFI0000AAAA", "cell by gene matrix")
    second = _file("IGVFFI0000BBBB", "cell by gene matrix")
    for accession in ("IGVFFI0000AAAA", "IGVFFI0000BBBB"):
        seed_record(
            connection,
            _file(accession, "cell by gene matrix", filtered=True),
            frame="object",
        )
    with pytest.raises(ValueError, match="none satisfied"):
        _content_type_record(connection, _analysis_set(first, second))


def test_write_download_entry(connection: PConnection) -> None:
    """Each entry is a URL then an indented out= option, as aria2c --input-file expects."""
    f_out = io.StringIO()
    get_url._write_download_entry(
        record=_file("IGVFFI0000FFFF", "fragments"),
        content_type=ContentType.FRAGMENTS,
        accession=AccessionId("IGVFDS0000AAAA"),
        connection=connection,
        f_out=f_out,
        logger=_LOGGER,
    )
    assert f_out.getvalue() == (
        "https://igvf-public.s3.us-west-2.amazonaws.com/2025/01/01/uuid/IGVFFI0000FFFF.h5ad\n"
        "  out=IGVFDS0000AAAA.bed.gz\n"
    )


def test_get_url_analysis_set(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, connection: PConnection
) -> None:
    """An "<input>;<output>" accession names the downloads after the output accession."""
    fragments = _file("IGVFFI0000FFFF", "fragments")
    analysis_set = _analysis_set(_file("IGVFFI0000MMMM", "cell by gene matrix"), fragments)
    seed_record(connection, analysis_set)
    seed_record(connection, fragments)
    monkeypatch.setattr(get_url.PConnection, "new", lambda *_a, **_k: connection)
    output = tmp_path / "aria2c_input.txt"
    get_url.get_url("IGVFDS0000AAAA;IGVFDS0000ZZZZ", output=output)
    # a second accession appends to the same input file
    get_url.get_url("IGVFFI0000FFFF", output=output)
    assert output.read_text().splitlines()[1::2] == [
        "  out=IGVFDS0000ZZZZ.h5ad",
        "  out=IGVFDS0000ZZZZ.bed.gz",
        # a file's downloads are named after its file set
        "  out=IGVFDS0000AAAA.bed.gz",
    ]


def test_get_url_file_wrong_format(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, connection: PConnection
) -> None:
    """A file accession given directly must be in a format the pipeline can use."""
    seed_record(connection, _file("IGVFFI0000MMMM", "cell by gene matrix", file_format="hdf5"))
    monkeypatch.setattr(get_url.PConnection, "new", lambda *_a, **_k: connection)
    with pytest.raises(ValueError, match="in format hdf5, but only h5ad"):
        get_url.get_url("IGVFFI0000MMMM", output=tmp_path / "out.txt")


def test_get_url_file_invalidated(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, connection: PConnection
) -> None:
    """A file named directly is downloaded even if invalidated: fix_bad_beds.nf needs them."""
    seed_record(
        connection,
        _file("IGVFFI0000BBBB", "peaks", upload_status="invalidated", file_set=_PSEUDOBULK_SET),
    )
    monkeypatch.setattr(get_url.PConnection, "new", lambda *_a, **_k: connection)
    output = tmp_path / "out.txt"
    get_url.get_url("IGVFFI0000BBBB;IGVFFI0000BBBB", output=output)
    assert output.read_text().splitlines()[1] == "  out=IGVFFI0000BBBB.tsv.gz"


def test_get_url_unsupported_content_type(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, connection: PConnection
) -> None:
    seed_record(connection, _file("IGVFFI0000QQQQ", "per-cell quality report"))
    monkeypatch.setattr(get_url.PConnection, "new", lambda *_a, **_k: connection)
    with pytest.raises(ValueError, match="per-cell quality report"):
        get_url.get_url("IGVFFI0000QQQQ", output=tmp_path / "out.txt")


def _reference(accession: str, summary: str) -> IgvfRecord:
    return cast(
        IgvfRecord,
        {
            "@id": f"/reference-files/{accession}/",
            "@type": ["ReferenceFile", "File", "Item"],
            "accession": accession,
            "aliases": [],
            "content_type": "genome reference",
            "status": "released",
            "file_set": {"summary": summary},
        },
    )


def test_get_species(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    connection: PConnection,
) -> None:
    """Species comes from the references of the analysis sets' matrix and fragments files."""
    seed_record(connection, _reference("IGVFFI0000GGGG", "Homo sapiens GRCh38 genome"))
    seed_record(connection, _reference("IGVFFI0000TTTT", "Homo sapiens GRCh38 transcriptome"))
    matrix = _file(
        "IGVFFI0000MMMM",
        "cell by gene matrix",
        reference_files=["/reference-files/IGVFFI0000TTTT/"],
    )
    fragments = _file(
        "IGVFFI0000FFFF",
        "fragments",
        reference_files=[{"@id": "/reference-files/IGVFFI0000GGGG/"}],
    )
    for record in (matrix, fragments, _analysis_set(matrix, fragments)):
        seed_record(connection, record)
    monkeypatch.setattr(get_species.PConnection, "new", lambda *_a, **_k: connection)
    get_species.get_species("IGVFDS0000AAAA, IGVFFI0000FFFF")
    assert capsys.readouterr().out == "human\n"


@pytest.mark.parametrize(
    ("summary", "species"),
    [
        ("Homo sapiens GRCh38 genome", "human"),
        ("Homo sapiens custom genome", "human"),
        ("Mus musculus GRCm39 mm10 GENCODE M25 transcriptome", "mouse"),
    ],
)
def test_get_summary_species(summary: str, species: str) -> None:
    assert get_species._get_summary_species(summary) == species


@pytest.mark.parametrize(
    "summary",
    [
        # a mixed human-mouse ("barnyard") reference names neither species
        "GRCh38, mm10 GENCODE 32, GENCODE M23 genome",
        "Homo sapiens and Mus musculus genome",
        "Danio rerio GRCz11 genome",
    ],
)
def test_get_summary_species_unknown(summary: str) -> None:
    with pytest.raises(ValueError, match="Unable to determine species"):
        get_species._get_summary_species(summary)


def test_get_species_mixed(monkeypatch: pytest.MonkeyPatch, connection: PConnection) -> None:
    seed_record(connection, _reference("IGVFFI0000HHHH", "Homo sapiens GRCh38 genome"))
    seed_record(connection, _reference("IGVFFI0000MOUS", "Mus musculus GRCm39 genome"))
    seed_record(
        connection,
        _file(
            "IGVFFI0000AAAA",
            "fragments",
            reference_files=["/reference-files/IGVFFI0000HHHH/"],
        ),
    )
    seed_record(
        connection,
        _file(
            "IGVFFI0000BBBB",
            "fragments",
            reference_files=["/reference-files/IGVFFI0000MOUS/"],
        ),
    )
    monkeypatch.setattr(get_species.PConnection, "new", lambda *_a, **_k: connection)
    with pytest.raises(ValueError, match="Found multiple species"):
        get_species.get_species("IGVFFI0000AAAA,IGVFFI0000BBBB")


@pytest.mark.parametrize(
    ("submitted_file_name", "expected"),
    [
        ("pseudobulks/x/peaks.narrowPeak.gz", "pseudobulks/x/peaks.narrowPeak.bb"),
        ("pseudobulks/x/fragments.tsv.gz", "pseudobulks/x/fragments.bb"),
        ("pseudobulks/x/peaks.bed.gz", "pseudobulks/x/peaks.bb"),
        ("pseudobulks/x/peaks.bed", "pseudobulks/x/peaks.bb"),
        ("peaks", "peaks.bb"),
    ],
)
def test_bed_record_to_bigbed_record_name(submitted_file_name: str, expected: str) -> None:
    bed_record = _file(
        "IGVFFI0000BBBB",
        "peaks",
        aliases=["lab:pseudobulk-x-peaks_narrowPeak_gz"],
        submitted_file_name=submitted_file_name,
    )
    assert find_bad_beds._bed_record_to_bigbed_record(bed_record)["submitted_file_name"] == expected


@pytest.mark.parametrize(
    ("content_type", "alias", "file_format_type", "bb_alias"),
    [
        ("peaks", "lab:x-peaks_bed_gz", "bed6+", "lab:x-peaks_bb"),
        ("fragments", "lab:x-fragments_tsv_gz", "bed3+", "lab:x-fragments_bb"),
        ("peaks", "lab:x-peaks", "bed6+", "lab:x-peaks_bb"),
    ],
)
def test_bed_record_to_bigbed_record(
    content_type: str, alias: str, file_format_type: str, bb_alias: str
) -> None:
    bed_record = _file(
        "IGVFFI0000BBBB",
        content_type,
        aliases=[alias],
        file_format="bed",
        submitted_file_name="x.bed.gz",
    )
    bb_record = find_bad_beds._bed_record_to_bigbed_record(bed_record)
    assert bb_record["file_format"] == "bigBed"
    assert bb_record["file_format_type"] == file_format_type
    assert bb_record["aliases"] == [bb_alias]
    # the BED record itself is left alone
    assert bed_record["file_format"] == "bed"
    assert bed_record["aliases"] == [alias]


def _bed_with_error(error: str) -> IgvfRecord:
    return _file(
        "IGVFFI0000BBBB",
        "peaks",
        aliases=["lab:x-peaks_bed_gz"],
        submitted_file_name="x/peaks.bed.gz",
        reference_files=["/reference-files/IGVFFI0000GGGG/"],
        validation_error_detail=json.dumps({"validate_files": error}),
    )


@pytest.fixture
def bad_beds_worker(tmp_path: Path, connection: PConnection) -> tuple[Path, PConnection]:
    """Set up find_bad_beds' thread-local state, with a genome reference to find."""
    seed_record(
        connection,
        _reference("IGVFFI0000GGGG", "Homo sapiens GRCh38 genome"),
        frame="object",
    )
    find_bad_beds._THREAD_LOCAL_DATA.output = tmp_path
    return tmp_path, connection


@pytest.mark.parametrize(
    "error",
    [
        "Error: bed->chromEnd[248956423] > chromSize[248956422]",
        "Error: score (1001) must be between 0 and 1000",
    ],
)
def test_check_bed_needs_filter(bad_beds_worker: tuple[Path, PConnection], error: str) -> None:
    output, connection = bad_beds_worker
    big_bed = _file("IGVFFI0000CCCC", "peaks", file_format="bigBed")
    find_bad_beds._check_bed(
        bed_record=_bed_with_error(error),
        big_bed_record=big_bed,
        content_type="peaks",
        pseudobulk_id=PortalId("/pseudobulk-sets/IGVFDS0000PPPP/"),
        connection=connection,
        logger=connection.logger,
    )
    assert (output / "IGVFFI0000BBBB.ref.txt").read_text() == "IGVFFI0000GGGG"
    assert json.loads((output / "IGVFFI0000BBBB.bed.json").read_text())["accession"] == (
        "IGVFFI0000BBBB"
    )
    # an existing big-bed record is re-used
    assert json.loads((output / "IGVFFI0000BBBB.bb.json").read_text())["accession"] == (
        "IGVFFI0000CCCC"
    )


def test_check_bed_ok(bad_beds_worker: tuple[Path, PConnection]) -> None:
    """A valid BED with a big-bed needs nothing."""
    output, connection = bad_beds_worker
    find_bad_beds._check_bed(
        bed_record=_bed_with_error("Error: chrom KI270728.1 not found"),
        big_bed_record=_file("IGVFFI0000CCCC", "peaks", file_format="bigBed"),
        content_type="peaks",
        pseudobulk_id=PortalId("/pseudobulk-sets/IGVFDS0000PPPP/"),
        connection=connection,
        logger=connection.logger,
    )
    assert list(output.iterdir()) == []


def test_check_bed_missing_big_bed(bad_beds_worker: tuple[Path, PConnection]) -> None:
    """A valid BED without a big-bed gets a new big-bed record."""
    output, connection = bad_beds_worker
    find_bad_beds._check_bed(
        bed_record=_bed_with_error(""),
        big_bed_record=None,
        content_type="peaks",
        pseudobulk_id=PortalId("/pseudobulk-sets/IGVFDS0000PPPP/"),
        connection=connection,
        logger=connection.logger,
    )
    bb_record = json.loads((output / "IGVFFI0000BBBB.bb.json").read_text())
    assert bb_record["file_format"] == "bigBed"
    assert bb_record["submitted_file_name"] == "x/peaks.bb"


def _payloads(infile: Path, drop_extra_fields: bool = False) -> list[dict[str, object]]:
    return list(
        register._iter_payloads(
            schema=cast(register.IgvfSchema, FAKE_SCHEMA),
            infile=infile,
            drop_extra_fields=drop_extra_fields,
        )
    )


def test_iter_payloads_from_tsv(tmp_path: Path) -> None:
    infile = tmp_path / "tabular_file.4.tsv"
    infile.write_text(
        "aliases\tfile_size\tcontrolled_access\tdescription\tattachment\tderived_from\t#note\n"
        'lab:a\t10\tFalse\tCD4 "naive"\t{"path": "x.pdf"}\t/files/A/,/files/B/\tignored\n'
        "lab:b\t\tTrue\t\t\t\t\n"
    )
    assert _payloads(infile) == [
        {
            "_profile": "tabular_file",
            "aliases": ["lab:a"],
            "file_size": 10,
            "controlled_access": False,
            "description": 'CD4 "naive"',
            "attachment": {"path": "x.pdf"},
            "derived_from": ["/files/A/", "/files/B/"],
        },
        # empty cells are left out of the payload
        {"_profile": "tabular_file", "aliases": ["lab:b"], "controlled_access": True},
    ]


def test_iter_payloads_from_tsv_extra_field(tmp_path: Path) -> None:
    infile = tmp_path / "x.tsv"
    infile.write_text("aliases\tunknown\nlab:a\t1\n")
    with pytest.raises(ValueError, match="Unknown field name 'unknown'"):
        _payloads(infile)
    assert _payloads(infile, drop_extra_fields=True) == [
        {"_profile": "tabular_file", "aliases": ["lab:a"]}
    ]


def test_iter_payloads_from_concatenated_json(tmp_path: Path) -> None:
    """UPLOAD_FIXED_BEDS concatenates one JSON record per file; embedded links become @ids."""
    infile = tmp_path / "tabular_files.json"
    first = {
        "aliases": ["lab:a"],
        "file_set": {"@id": "/pseudobulk-sets/X/", "summary": "s"},
    }
    second = {"aliases": ["lab:b"], "derived_from": [{"@id": "/files/A/"}], "uuid": "u"}
    infile.write_text(json.dumps(first) + json.dumps(second) + "\n")
    assert _payloads(infile, drop_extra_fields=True) == [
        {
            "aliases": ["lab:a"],
            "file_set": "/pseudobulk-sets/X/",
            "_profile": "tabular_file",
        },
        {
            "aliases": ["lab:b"],
            "derived_from": ["/files/A/"],
            "_profile": "tabular_file",
        },
    ]
    with pytest.raises(ValueError, match="uuid"):
        _payloads(infile, drop_extra_fields=False)


def test_iter_payloads_unknown_format(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="Unknown payload specification format"):
        _payloads(tmp_path / "x.csv")


class _PatchConnection:
    IGVFID_KEY = "@id"

    def __init__(self) -> None:
        self.patched: list[dict[str, object]] = []

    # extend_array_values is passed by keyword, so it must keep its name
    def patch(self, payload: dict[str, object], extend_array_values: bool) -> None:  # noqa: ARG002
        self.patched.append(payload)


def test_patch_payload_without_aliases() -> None:
    """A patch is identified by record_id, so a payload without aliases is fine."""
    fake_connection = _PatchConnection()
    register._THREAD_LOCAL_DATA.connection = fake_connection
    register._THREAD_LOCAL_DATA.register_config = register.RegisterConfig(
        igvf_mode=IgvfMode.prod, profile_id="tabular_file"
    )
    record_id = register._patch_payload({"record_id": "IGVFFI0000AAAA", "description": "d"})
    assert record_id == "IGVFFI0000AAAA"
    assert fake_connection.patched == [{"description": "d", "@id": "IGVFFI0000AAAA"}]


def test_patch_payload_requires_record_id() -> None:
    register._THREAD_LOCAL_DATA.connection = _PatchConnection()
    with pytest.raises(ValueError, match="record_id"):
        register._patch_payload({"description": "d"})
