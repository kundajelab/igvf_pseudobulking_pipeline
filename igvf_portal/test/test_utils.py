import gzip
import logging
from pathlib import Path
from typing import cast

import pytest
import requests

from igvf_portal import utils
from igvf_portal.types import GeneInfoRow, IgvfRecord


@pytest.mark.parametrize("name", ["table.tsv", "table.tsv.gz"])
def test_maybe_gzipped_round_trip(tmp_path: Path, name: str) -> None:
    """Text written through maybe_gzipped is all there when read back, compressed or not."""
    path = tmp_path / name
    text = "a\tb\n" * 1000
    with utils.maybe_gzipped(path, mode="w") as f_out:
        f_out.write(text)
    with utils.maybe_gzipped(path, mode="r") as f_in:
        assert f_in.read() == text


def test_maybe_gzipped_writes_valid_gzip(tmp_path: Path) -> None:
    path = tmp_path / "table.tsv.gz"
    with utils.maybe_gzipped(path, mode="w") as f_out:
        f_out.write("x\n")
    assert gzip.decompress(path.read_bytes()) == b"x\n"


@pytest.mark.parametrize("name", ["gene_info.csv", "gene_info.csv.gz"])
def test_write_csv(tmp_path: Path, name: str) -> None:
    rows: list[GeneInfoRow] = [
        {"gene_id": "ENSG1", "gene_name": "MT-ND1", "mt": True, "ribo": False},
        {"gene_id": "ENSG2", "gene_name": "RPL1", "mt": False, "ribo": True},
    ]
    output = tmp_path / name
    utils.write_csv(rows=rows, output_csv=output, row_type=GeneInfoRow)
    assert list(utils.iter_csv_rows(output)) == [
        {"gene_id": "ENSG1", "gene_name": "MT-ND1", "mt": "True", "ribo": "False"},
        {"gene_id": "ENSG2", "gene_name": "RPL1", "mt": "False", "ribo": "True"},
    ]


def test_iter_csv_rows_missing_columns(tmp_path: Path) -> None:
    path = tmp_path / "x.tsv"
    path.write_text("a\tb\n1\t2\n")
    with pytest.raises(ValueError, match="missing required columns: c"):
        list(utils.iter_csv_rows(path, required_columns=("a", "c")))


@pytest.mark.parametrize(
    "file_name",
    [
        "/repo/igvf_portal/src/igvf_portal/tools/get_url.py",
        "/opt/env/lib/python3.14t/site-packages/igvf_portal/tools/get_url.py",
        "/home/user/.venv/igvf_portal/tools/get_url.py",
    ],
)
def test_get_logger_from_file(file_name: str) -> None:
    """The logger name comes from the package and tool, wherever the package is installed."""
    assert utils.get_logger_from_file(file_name).name == "igvf-portal get-url"


def test_get_sep_from_path() -> None:
    assert utils.get_sep_from_path(Path("a.tsv.gz")) == "\t"
    assert utils.get_sep_from_path(Path("a.b.csv")) == ","
    with pytest.raises(ValueError, match="Could not infer separator"):
        utils.get_sep_from_path(Path("a.txt"))


def test_parse_s3_uri() -> None:
    assert utils.parse_s3_uri("s3://bucket/a/b/c.txt") == ("bucket", "a/b/c.txt")
    with pytest.raises(ValueError, match="not a valid S3 path"):
        utils.parse_s3_uri("https://bucket/a")


def test_sanitize_to_ascii_underscore() -> None:
    assert utils.sanitize_to_ascii_underscore("Café, a.b/c\\d  e") == "Cafe_a_b_c_d_e"


class _Sleeps:
    def __init__(self) -> None:
        self.calls: list[float] = []

    def __call__(self, seconds: float) -> None:
        self.calls.append(seconds)


@pytest.fixture
def sleeps(monkeypatch: pytest.MonkeyPatch) -> _Sleeps:
    """Record retry delays instead of sleeping."""
    recorder = _Sleeps()
    monkeypatch.setattr(utils.time, "sleep", recorder)
    return recorder


def test_retry_succeeds_after_failures(sleeps: _Sleeps) -> None:
    attempts: list[int] = []

    def flaky() -> str:
        attempts.append(1)
        if len(attempts) < 3:
            raise ConnectionError("try again")
        return "ok"

    assert utils.retry(num_tries=3, delay=1.0, backoff=2.0)(flaky)() == "ok"
    assert len(attempts) == 3
    assert sleeps.calls == [1.0, 2.0]


def test_retry_delay_restarts_each_call(sleeps: _Sleeps) -> None:
    """The backoff of one call does not carry into the next call of the same function."""

    def always_fails() -> None:
        raise ConnectionError("no")

    retried = utils.retry(num_tries=2, delay=1.0, backoff=2.0)(always_fails)
    for _ in range(3):
        with pytest.raises(ConnectionError):
            retried()
    assert sleeps.calls == [1.0, 1.0, 1.0]


def test_retry_does_not_retry_subclass_of_no_retry_exception(
    sleeps: _Sleeps, caplog: pytest.LogCaptureFixture
) -> None:
    """Requests' JSONDecodeError subclasses json's, and must not be retried either."""
    import json

    def bad_json() -> None:
        raise requests.exceptions.JSONDecodeError("bad", "doc", 0)

    retried = utils.retry(
        num_tries=3,
        no_retry_exceptions=(json.JSONDecodeError,),
        logger=logging.getLogger("test"),
    )(bad_json)
    with pytest.raises(requests.exceptions.JSONDecodeError):
        retried()
    assert sleeps.calls == []
    assert "Not retrying exception of type 'JSONDecodeError'" in caplog.text


def test_retry_rejects_no_tries() -> None:
    with pytest.raises(ValueError, match="num_tries"):
        utils.retry(num_tries=0)


def test_file_not_in_portal() -> None:
    assert utils.file_not_in_portal(cast(IgvfRecord, {"status": "deleted"}))
    assert utils.file_not_in_portal(
        cast(IgvfRecord, {"status": "in progress", "upload_status": "file not found"})
    )
    assert not utils.file_not_in_portal(
        cast(IgvfRecord, {"status": "released", "upload_status": "validated"})
    )


def test_retry_should_retry(sleeps: _Sleeps) -> None:
    """should_retry can tell errors of one type apart, e.g. a temporary from a permanent one."""
    attempts: list[str] = []

    def fails(message: str) -> None:
        attempts.append(message)
        raise ValueError(message)

    retried = utils.retry(num_tries=3, should_retry=lambda e: "temporary" in str(e))(fails)
    with pytest.raises(ValueError, match="permanent"):
        retried("permanent")
    assert attempts == ["permanent"]
    with pytest.raises(ValueError, match="temporary"):
        retried("temporary")
    assert attempts == ["permanent"] + ["temporary"] * 3
    assert sleeps.calls == [5.0, 10.0]
