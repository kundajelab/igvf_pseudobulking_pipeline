import gzip
import json
import urllib.parse
from collections.abc import Callable, Collection, Iterator
from pathlib import Path
from typing import Final, cast

import pytest
import requests
from igvf_utils.profiles import IgvfSchema
from requests.adapters import HTTPAdapter

from igvf_portal import connection as connection_module
from igvf_portal import utils
from igvf_portal.connection import PConnection, _local_accession, _query_values
from igvf_portal.enums import IgvfMode
from igvf_portal.register_config import RegisterConfig
from igvf_portal.tools import register
from igvf_portal.types import AccessionId, Alias, IgvfRecord, PortalId, RecordNotFound

from .conftest import seed_record


@pytest.mark.parametrize("igvf_mode", list(IgvfMode))
def test_inherited_methods_use_same_portal(igvf_mode: IgvfMode) -> None:
    """Inherited igvf_utils methods and PConnection's own methods must use the same Portal.

    Methods inherited from igvf_utils (patch, profiles, upload credentials) use igvf_mode.url, and
    those defined in PConnection use mode.url. They must name the same Portal: otherwise a run could
    POST to one Portal but PATCH and fetch upload credentials from another.
    """
    # NOTE: no request is made: igvf_mode is not validated (and the network is blocked in tests)
    connection = PConnection.new(igvf_mode)
    inherited_host = urllib.parse.urlparse(connection.igvf_mode.url).netloc
    assert inherited_host == urllib.parse.urlparse(connection.mode.url).netloc


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("IGVFFI0000AAAA", "IGVFFI0000AAAA"),
        ("/tabular-files/IGVFFI0000AAAA/", "IGVFFI0000AAAA"),
        ("/IGVFSM0104SYAP/", "IGVFSM0104SYAP"),
        ("anshul-kundaje:some-alias", None),
        ("/sample-terms/CL_0000084/", None),
        ("IGVFFI0000AAAAB", None),
    ],
)
def test_local_accession(value: str, expected: str | None) -> None:
    assert _local_accession(value) == expected


def test_query_values() -> None:
    assert _query_values(None) == []
    assert _query_values("a") == ["a"]
    assert _query_values(["a", "b"]) == ["a", "b"]
    assert _query_values(True) == ["true"]
    assert _query_values([False, 3]) == ["false", "3"]


def test_search_url(connection: PConnection) -> None:
    url = connection._search_url(
        "PseudobulkSet",
        field_filters={"status!": "deleted", "lab.@id": ["/labs/a/", "/labs/b/"]},
        field=("@id",),
        limit=5,
    )
    parsed = urllib.parse.urlparse(url)
    assert f"{parsed.scheme}://{parsed.netloc}{parsed.path}" == (
        "https://api.data.igvf.org/search/"
    )
    assert urllib.parse.parse_qsl(parsed.query) == [
        ("format", "json"),
        ("limit", "5"),
        ("type", "PseudobulkSet"),
        ("field", "@id"),
        ("status!", "deleted"),
        ("lab.@id", "/labs/a/"),
        ("lab.@id", "/labs/b/"),
        ("frame", "object"),
    ]


def test_search_records_rejects_long_url(connection: PConnection) -> None:
    with pytest.raises(ValueError, match="search_records_matching_any"):
        connection.search_records(
            "File",
            field_filters={"accession": [f"IGVFFI{i:04d}AAAA" for i in range(500)]},
        )


def test_search_records_matching_any_batches(
    monkeypatch: pytest.MonkeyPatch, connection: PConnection
) -> None:
    """A long list of values is split into searches that fit, and every value is searched once."""
    searched: list[list[str]] = []

    def fake_search_records(
        _type: str, *, field_filters: dict[str, object], **_kwargs: object
    ) -> list[dict[str, object]]:
        batch = cast(list[str], field_filters["accession"])
        url = connection._search_url(_type, field_filters=field_filters)
        assert len(url) <= connection_module._MAX_SEARCH_URL_LENGTH
        searched.append(batch)
        # the same record is found by more than one batch, and is returned only once
        return [{"@id": f"/files/{value}/"} for value in batch] + [{"@id": "/files/shared/"}]

    monkeypatch.setattr(connection, "search_records", fake_search_records)
    values = [f"IGVFFI{i:04d}AAAA" for i in range(500)]
    records = connection.search_records_matching_any("File", "accession", values + values[:10])
    assert len(searched) > 1
    assert sorted(value for batch in searched for value in batch) == sorted(values)
    assert len(records) == len(values) + 1


def test_search_records_matching_any_no_values(connection: PConnection) -> None:
    assert connection.search_records_matching_any("File", "accession", []) == []


def _record(accession: str, **fields: object) -> IgvfRecord:
    return cast(
        IgvfRecord,
        {
            "@id": f"/tabular-files/{accession}/",
            "accession": accession,
            "aliases": [f"lab:{accession.lower()}"],
            "status": "released",
        }
        | fields,
    )


def test_lookup_record_uses_cache_for_every_identifier(connection: PConnection) -> None:
    """A cached record is found by its @id, accession, or alias, for the frame it was cached in."""
    record = _record("IGVFFI0000AAAA")
    seed_record(connection, record, frame="object")
    for key in (
        "/tabular-files/IGVFFI0000AAAA/",
        "IGVFFI0000AAAA",
        "lab:igvffi0000aaaa",
    ):
        assert connection.lookup_record(Alias(key), frame="object") is record


def test_lookup_record_caches_get(monkeypatch: pytest.MonkeyPatch, connection: PConnection) -> None:
    calls: list[object] = []

    def fake_get(rec_ids: list[str], **_kwargs: object) -> IgvfRecord:
        calls.append(rec_ids)
        return _record("IGVFFI0000AAAA")

    monkeypatch.setattr(connection, "get", fake_get)
    connection.lookup_record(AccessionId("IGVFFI0000AAAA"))
    connection.lookup_record(PortalId("/tabular-files/IGVFFI0000AAAA/"))
    assert len(calls) == 1


def test_lookup_record_negative_cache(
    monkeypatch: pytest.MonkeyPatch, connection: PConnection
) -> None:
    calls: list[object] = []

    def fake_get(rec_ids: list[str], **_kwargs: object) -> None:
        calls.append(rec_ids)

    monkeypatch.setattr(connection, "get", fake_get)
    for _ in range(2):
        with pytest.raises(RecordNotFound):
            connection.lookup_record(Alias("lab:missing"), cache_negative=True)
    assert len(calls) == 1
    # without cache_negative the record is asked for again
    with pytest.raises(RecordNotFound):
        connection.lookup_record(Alias("lab:other"))
    with pytest.raises(RecordNotFound):
        connection.lookup_record(Alias("lab:other"))
    assert len(calls) == 3


class _Response:
    def __init__(self, status_code: int, body: object = None) -> None:
        self.status_code = status_code
        self.ok = status_code < 400
        self.body = body

    def json(self) -> object:
        return self.body

    def raise_for_status(self) -> None:
        if not self.ok:
            raise requests.HTTPError(f"{self.status_code}", response=None)


def _fake_session_get(responses: dict[str, _Response]) -> Callable[..., _Response]:
    def fake_get(url: str, **_kwargs: object) -> _Response:
        for key, response in responses.items():
            if f"/{key}/" in url or url.endswith(f"/{key}?format=json"):
                return response
        raise AssertionError(f"Unexpected URL {url}")

    return fake_get


def test_get_tries_each_identifier(
    monkeypatch: pytest.MonkeyPatch, connection: PConnection
) -> None:
    record = {"@id": "/files/IGVFFI0000AAAA/"}
    monkeypatch.setattr(
        connection.session,
        "get",
        _fake_session_get({"lab:a": _Response(404), "IGVFFI0000AAAA": _Response(200, record)}),
    )
    assert connection.get(["lab:a", "IGVFFI0000AAAA"]) == record


def test_get_not_found(monkeypatch: pytest.MonkeyPatch, connection: PConnection) -> None:
    monkeypatch.setattr(connection.session, "get", _fake_session_get({"lab:a": _Response(404)}))
    assert connection.get("lab:a", ignore404=True) is None
    with pytest.raises(RecordNotFound):
        connection.get("lab:a", ignore404=False)


def test_get_forbidden(monkeypatch: pytest.MonkeyPatch, connection: PConnection) -> None:
    monkeypatch.setattr(connection.session, "get", _fake_session_get({"lab:a": _Response(403)}))
    with pytest.raises(RuntimeError, match="forbidden"):
        connection.get("lab:a")


def test_infer_principal_accessions(connection: PConnection) -> None:
    """Intermediate analysis sets lead to the principal sets they are input for."""
    principal = _record(
        "IGVFDS0000PPPP",
        input_file_sets=[
            {"accession": "IGVFDS0000IIII"},
            {"accession": "IGVFDS0000JJJJ"},
        ],
    )
    seed_record(connection, principal)
    seed_record(connection, _record("IGVFDS0000IIII", input_for=[principal["@id"]]))
    seed_record(connection, _record("IGVFDS0000JJJJ", input_for=[principal["@id"]]))
    # a set that is not input for anything is itself principal; the Portal omits empty input_for
    seed_record(connection, _record("IGVFDS0000KKKK"))
    seed_record(connection, _record("IGVFDS0000DDDD", status="deleted"))
    principal_ids = connection.infer_principal_accessions(
        [
            AccessionId("IGVFDS0000IIII"),
            AccessionId("IGVFDS0000JJJJ"),
            AccessionId("IGVFDS0000KKKK"),
            AccessionId("IGVFDS0000DDDD"),
        ]
    )
    assert principal_ids == {principal["@id"], "/tabular-files/IGVFDS0000KKKK/"}


@pytest.mark.parametrize(
    ("key", "existing", "new", "equal"),
    [
        ("description", "a", "a", True),
        ("description", "a", "b", False),
        ("aliases", ["lab:b", "lab:a"], "lab:a,lab:b", True),
        ("file_set", "/analysis-sets/IGVFDS0000AAAA/", "IGVFDS0000AAAA", True),
        ("derived_from", ["/files/IGVFFI0000AAAA/"], "IGVFFI0000AAAA", True),
        ("derived_from", ["/files/IGVFFI0000AAAA/"], "IGVFFI0000BBBB", False),
        ("file_size", 10, 11, False),
    ],
)
def test_props_equal(
    connection: PConnection, key: str, existing: object, new: object, equal: bool
) -> None:
    """Compare values offline: accessions and @ids ending in one need no Portal lookup."""
    assert connection._props_equal(key, existing, new) is equal


class _CredentialsError:
    def __init__(self, status_code: int) -> None:
        self.status_code = status_code

    def __call__(self, file_id: str) -> None:
        response = requests.Response()
        response.status_code = self.status_code
        raise requests.HTTPError(f"{self.status_code} for {file_id}", response=response)


def test_upload_file_skips_finalized(
    monkeypatch: pytest.MonkeyPatch, connection: PConnection, tmp_path: Path
) -> None:
    """The Portal answers 403 when a file is finalized, so there is nothing to upload."""
    monkeypatch.setattr(connection, "_regenerate_s3_client", _CredentialsError(403))
    connection.upload_file(AccessionId("IGVFFI0000AAAA"), file_path=tmp_path / "x.bed.gz")


@pytest.mark.parametrize("status_code", [401, 404, 500, 503])
def test_upload_file_raises_on_other_credential_errors(
    monkeypatch: pytest.MonkeyPatch, connection: PConnection, tmp_path: Path, status_code: int
) -> None:
    """Any other failure to get credentials leaves the file un-uploaded, so it must not pass."""
    monkeypatch.setattr(connection, "_regenerate_s3_client", _CredentialsError(status_code))
    with pytest.raises(requests.HTTPError, match=str(status_code)):
        connection.upload_file(AccessionId("IGVFFI0000AAAA"), file_path=tmp_path / "x.bed.gz")


def test_upload_file_raises_on_forbidden_without_continue(
    monkeypatch: pytest.MonkeyPatch, connection: PConnection, tmp_path: Path
) -> None:
    connection.continue_on_failed_credentials = False
    monkeypatch.setattr(connection, "_regenerate_s3_client", _CredentialsError(403))
    with pytest.raises(requests.HTTPError, match="403"):
        connection.upload_file(AccessionId("IGVFFI0000AAAA"), file_path=tmp_path / "x.bed.gz")


class _FakeResponse:
    """Just enough of requests.Response for the download methods."""

    def __init__(
        self,
        url: str,
        status_code: int = 200,
        history: list[_FakeResponse] | None = None,
        chunks: list[bytes] | None = None,
    ) -> None:
        self.url = url
        self.status_code = status_code
        self.history = history or []
        self.chunks = chunks or []

    @property
    def content(self) -> bytes:
        return b"".join(self.chunks)

    def __enter__(self) -> _FakeResponse:
        return self

    def __exit__(self, *_args: object) -> None:
        pass

    def raise_for_status(self) -> None:
        if self.status_code >= 400:
            raise requests.HTTPError(f"{self.status_code}")

    def iter_content(self, chunk_size: int) -> Iterator[bytes]:
        data = self.content
        for start in range(0, len(data), chunk_size):
            yield data[start : start + chunk_size]


def _download_record(**fields: str) -> IgvfRecord:
    return cast(
        IgvfRecord,
        {
            "accession": "IGVFFI0000AAAA",
            "href": "/tabular-files/IGVFFI0000AAAA/@@download/x",
        }
        | fields,
    )


def test_get_record_http_url_public(connection: PConnection) -> None:
    record = _download_record(s3_uri="s3://igvf-public/2025/01/01/uuid/IGVFFI0000AAAA.bed.gz")
    assert connection.get_record_http_url(record) == (
        "https://igvf-public.s3.us-west-2.amazonaws.com/2025/01/01/uuid/IGVFFI0000AAAA.bed.gz"
    )


def test_get_record_http_url_private_redirect(
    monkeypatch: pytest.MonkeyPatch, connection: PConnection
) -> None:
    """The presigned URL is returned even though S3 refuses the HEAD (it is signed for GET)."""
    presigned = "https://igvf-private.s3.amazonaws.com/x?X-Amz-Signature=abc"

    def fake_head(url: str, **kwargs: object) -> _FakeResponse:
        assert url == "https://api.data.igvf.org/tabular-files/IGVFFI0000AAAA/@@download/x"
        assert kwargs["timeout"] is not None
        return _FakeResponse(presigned, 403, history=[_FakeResponse(url, 307)])

    monkeypatch.setattr(connection.session, "head", fake_head)
    record = _download_record(s3_uri="s3://igvf-private/2025/x")
    assert connection.get_record_http_url(record) == presigned


@pytest.mark.parametrize("status_code", [200, 403, 404])
def test_get_record_http_url_private_not_redirected(
    monkeypatch: pytest.MonkeyPatch, connection: PConnection, status_code: int
) -> None:
    """If the Portal does not redirect, there is no download URL to hand to aria2c."""
    monkeypatch.setattr(
        connection.session, "head", lambda url, **_kwargs: _FakeResponse(url, status_code)
    )
    with pytest.raises((requests.HTTPError, RuntimeError)):
        connection.get_record_http_url(_download_record())


def _serve(
    monkeypatch: pytest.MonkeyPatch, connection: PConnection, data: bytes, status_code: int = 200
) -> list[str]:
    """Answer the connection's GETs with data, and return the list of URLs requested."""
    urls: list[str] = []

    def fake_get(url: str, **_kwargs: object) -> _FakeResponse:
        urls.append(url)
        return _FakeResponse(url, status_code, chunks=[data])

    monkeypatch.setattr(connection.session, "get", fake_get)
    return urls


@pytest.mark.parametrize("make_dir", [True, False])
def test_download_record_to_folder(
    monkeypatch: pytest.MonkeyPatch, connection: PConnection, tmp_path: Path, make_dir: bool
) -> None:
    """An output without a file suffix is a folder to download into, using the href name."""
    urls = _serve(monkeypatch, connection, b"abcdef")
    folder = tmp_path / "downloads"
    if make_dir:
        folder.mkdir()
    record = _download_record(href="/reference-files/IGVFFI0000AAAA/@@download/genome.fa.gz")
    output = connection.download_record(record, chunk_size=4, output=folder)
    assert output == folder / "genome.fa.gz"
    assert output.read_bytes() == b"abcdef"
    assert urls == ["https://api.data.igvf.org/IGVFFI0000AAAA/@@download"]


def test_download_record_to_file(
    monkeypatch: pytest.MonkeyPatch, connection: PConnection, tmp_path: Path
) -> None:
    _serve(monkeypatch, connection, b"abcdef")
    target = tmp_path / "sub" / "metadata.tsv"
    record = _download_record(href="/tabular-files/IGVFFI0000AAAA/@@download/other.tsv")
    assert connection.download_record(record, output=target) == target
    assert target.read_bytes() == b"abcdef"


@pytest.mark.parametrize("decompress", [True, False])
def test_stream_bytes_whole(
    monkeypatch: pytest.MonkeyPatch, connection: PConnection, decompress: bool
) -> None:
    """With chunk_size=None the content comes back as bytes, not an iterator."""
    compressed = gzip.compress(b"a\tb\n")
    _serve(monkeypatch, connection, compressed)
    content = connection.stream_bytes(
        AccessionId("IGVFFI0000AAAA"), decompress=decompress, chunk_size=None
    )
    assert content == (b"a\tb\n" if decompress else compressed)


@pytest.mark.parametrize("chunk_size", [1, 7, 100_000])
def test_stream_bytes_chunks_decompressed(
    monkeypatch: pytest.MonkeyPatch, connection: PConnection, chunk_size: int
) -> None:
    """Chunks split the gzip data anywhere, and bgzip writes several gzip members end to end."""
    text = b"".join(f"chr1\t{i}\t{i + 1}\n".encode() for i in range(2000))
    multi_member = gzip.compress(text[:7000]) + gzip.compress(text[7000:])
    _serve(monkeypatch, connection, multi_member)
    chunks = connection.stream_bytes(
        AccessionId("IGVFFI0000AAAA"), decompress=True, chunk_size=chunk_size
    )
    assert b"".join(chunks) == text


def test_stream_bytes_chunks_raw(monkeypatch: pytest.MonkeyPatch, connection: PConnection) -> None:
    _serve(monkeypatch, connection, b"abcdef")
    chunks = connection.stream_bytes(AccessionId("IGVFFI0000AAAA"), decompress=False, chunk_size=4)
    assert list(chunks) == [b"abcd", b"ef"]


def test_read_remote_bytes(monkeypatch: pytest.MonkeyPatch, connection: PConnection) -> None:
    _serve(monkeypatch, connection, gzip.compress(b"a\tb\n"))
    with connection.read_remote_bytes(AccessionId("IGVFFI0000AAAA")) as f_in:
        assert f_in.read() == "a\tb\n"


def test_read_remote_bytes_does_not_mask_caller_errors(
    monkeypatch: pytest.MonkeyPatch, connection: PConnection
) -> None:
    """An error in the with-block is the caller's, not a failure to read the bytes."""
    _serve(monkeypatch, connection, gzip.compress(b"a\tb\n"))
    with (
        pytest.raises(KeyError, match="caller"),
        connection.read_remote_bytes(AccessionId("IGVFFI0000AAAA")),
    ):
        raise KeyError("caller")


@pytest.mark.parametrize(("data", "status_code"), [(b"not gzip", 200), (b"", 404)])
def test_read_remote_bytes_reports_bad_download(
    monkeypatch: pytest.MonkeyPatch, connection: PConnection, data: bytes, status_code: int
) -> None:
    _serve(monkeypatch, connection, data, status_code=status_code)
    with (
        pytest.raises(RuntimeError, match="Could not open IGVFFI0000AAAA"),
        connection.read_remote_bytes(AccessionId("IGVFFI0000AAAA")),
    ):
        pass


def _http_error(status_code: int) -> requests.HTTPError:
    response = requests.Response()
    response.status_code = status_code
    return requests.HTTPError(f"{status_code}", response=response)


@pytest.mark.parametrize(
    ("exception", "transient"),
    [
        (_http_error(500), True),
        (_http_error(503), True),
        (_http_error(429), True),
        (_http_error(400), False),
        (_http_error(404), False),
        (requests.ConnectionError("reset"), True),
        (requests.exceptions.JSONDecodeError("bad", "<html>", 0), True),
        (RuntimeError("S3 upload failed"), True),
        (FileNotFoundError("x.bed.gz"), False),
    ],
)
def test_upload_error_is_transient(exception: Exception, transient: bool) -> None:
    assert connection_module._upload_error_is_transient(exception) is transient


class _FakeProfiles:
    """Just enough of igvf_utils' Profiles for guard_upload."""

    FILE_PROFILE_ID: Final[tuple[str, ...]] = ("tabular_file",)
    SUBMITTED_FILE_PROP_NAME: Final[str] = "submitted_file_name"


class _FlakyUpload:
    """Stand-in for upload_file that fails with the given errors, then succeeds."""

    def __init__(self, *errors: Exception) -> None:
        self.errors = list(errors)
        self.calls = 0

    def __call__(self, **_kwargs: object) -> None:
        self.calls += 1
        if self.errors:
            raise self.errors.pop(0)


def _guard_upload(
    monkeypatch: pytest.MonkeyPatch, connection: PConnection, upload: _FlakyUpload
) -> None:
    monkeypatch.setattr(connection, "_profiles", _FakeProfiles())
    monkeypatch.setattr(connection, "upload_file", upload)
    monkeypatch.setattr(connection_module.utils.time, "sleep", lambda _seconds: None)
    connection.guard_upload(
        upload_file=True,
        profile=cast(IgvfSchema, _FakeSchemaName("tabular_file")),
        accession_id=AccessionId("IGVFFI0000AAAA"),
        payload={"submitted_file_name": "x.bed.gz", "md5sum": "abc"},
    )


class _FakeSchemaName:
    def __init__(self, name: str) -> None:
        self.name = name


def test_guard_upload_retries_transient_errors(
    monkeypatch: pytest.MonkeyPatch, connection: PConnection
) -> None:
    upload = _FlakyUpload(_http_error(502), requests.ConnectionError("reset"))
    _guard_upload(monkeypatch, connection, upload)
    assert upload.calls == 3


@pytest.mark.parametrize("error", [_http_error(404), FileNotFoundError("x.bed.gz")])
def test_guard_upload_does_not_retry_permanent_errors(
    monkeypatch: pytest.MonkeyPatch, connection: PConnection, error: Exception
) -> None:
    upload = _FlakyUpload(error)
    with pytest.raises(type(error)):
        _guard_upload(monkeypatch, connection, upload)
    assert upload.calls == 1


def test_guard_upload_gives_up(monkeypatch: pytest.MonkeyPatch, connection: PConnection) -> None:
    upload = _FlakyUpload(*(_http_error(500) for _ in range(5)))
    with pytest.raises(requests.HTTPError, match="500"):
        _guard_upload(monkeypatch, connection, upload)
    assert upload.calls == connection.upload_retry.num_tries


def test_register_config_sets_upload_retry() -> None:
    config = RegisterConfig(
        igvf_mode=IgvfMode.prod, profile_id="tabular_file", num_tries=5, delay=1.0, backoff=3.0
    )
    assert config.new_connection.upload_retry == utils.RetryPolicy(
        num_tries=5, delay=1.0, backoff=3.0
    )


def test_session_retries(connection: PConnection) -> None:
    """The session retries every idempotent method the connection uses, plus POST."""
    adapter = connection.session.get_adapter("https://api.data.igvf.org")
    assert isinstance(adapter, HTTPAdapter)
    retry = adapter.max_retries
    assert set(cast(Collection[str], retry.allowed_methods)) == {"HEAD", "GET", "POST", "PUT"}
    assert retry.status == 3
    assert 502 in cast(Collection[int], retry.status_forcelist)


def _http_response(
    status_code: int, body: object = None, text: str | None = None, method: str = "POST"
) -> requests.Response:
    """A real requests.Response with the given JSON body, or else text (e.g. a gateway's HTML)."""
    response = requests.Response()
    response.status_code = status_code
    response.reason = "Reason"
    response.url = "https://api.data.igvf.org/x/"
    response.encoding = "utf-8"
    response._content = (json.dumps(body) if text is None else text).encode()
    response.request = requests.Request(method, response.url).prepare()
    return response


_GATEWAY_HTML: Final[str] = "<html><body><h1>502 Bad Gateway</h1></body></html>"


def test_response_detail() -> None:
    assert '"title": "Conflict"' in connection_module._response_detail(
        _http_response(409, {"title": "Conflict"})
    )
    assert connection_module._response_detail(_http_response(502, text=_GATEWAY_HTML)) == (
        f"HTTP 502: {_GATEWAY_HTML}"
    )


class _FakeProfile:
    name = "tabular_file"


def _prepare_post(
    monkeypatch: pytest.MonkeyPatch, connection: PConnection, response: requests.Response
) -> None:
    """Skip payload validation (which needs the Portal's schemas), and answer the POST."""
    payload: dict[str, object] = {"aliases": ["lab:a"]}
    monkeypatch.setattr(
        connection,
        "_prep_and_validate_payload",
        lambda **_kwargs: (payload, _FakeProfile(), False, ["lab:a"]),
    )
    monkeypatch.setattr(connection.session, "post", lambda *_args, **_kwargs: response)
    # posting logs to igvf_utils' post logger, which needs no network
    monkeypatch.setattr(connection, "_log_post", lambda **_kwargs: None)


def test_post_success(monkeypatch: pytest.MonkeyPatch, connection: PConnection) -> None:
    record = {"accession": "IGVFFI0000AAAA", "uuid": "u"}
    _prepare_post(monkeypatch, connection, _http_response(201, {"@graph": [record]}))
    assert connection.post({"aliases": ["lab:a"]}, upload_file=False) == record


@pytest.mark.parametrize("status_code", [500, 502, 504])
def test_post_gateway_error_is_http_error(
    monkeypatch: pytest.MonkeyPatch, connection: PConnection, status_code: int
) -> None:
    """A non-JSON error page raises HTTPError with its status, not a JSONDecodeError."""
    _prepare_post(monkeypatch, connection, _http_response(status_code, text=_GATEWAY_HTML))
    with pytest.raises(requests.HTTPError, match=str(status_code)):
        connection.post({"aliases": ["lab:a"]}, upload_file=False)


def test_patch_in_post_gateway_error_is_http_error(
    monkeypatch: pytest.MonkeyPatch, connection: PConnection
) -> None:
    monkeypatch.setattr(
        connection,
        "_get_record_changes",
        lambda **_kwargs: ({}, True, {}, {"description"}, set()),
    )
    monkeypatch.setattr(
        connection.session,
        "put",
        lambda *_args, **_kwargs: _http_response(502, text=_GATEWAY_HTML, method="PUT"),
    )
    with pytest.raises(requests.HTTPError, match="502"):
        connection._patch_in_post(
            record_id="/files/IGVFFI0000AAAA/",
            profile=cast(IgvfSchema, _FakeProfile()),
            payload={"description": "d"},
            existing_record={"description": "d"},
        )


def test_regenerate_aws_upload_creds(
    monkeypatch: pytest.MonkeyPatch, connection: PConnection
) -> None:
    """Credentials are requested through the session, so the request gets its retries."""
    credentials = {"access_key": "a", "secret_key": "s", "session_token": "t"}
    urls: list[str] = []

    def fake_post(url: str, **_kwargs: object) -> requests.Response:
        urls.append(url)
        return _http_response(200, {"@graph": [{"upload_credentials": credentials}]})

    monkeypatch.setattr(connection.session, "post", fake_post)
    assert connection.regenerate_aws_upload_creds("IGVFFI0000AAAA") == credentials
    assert urls == ["https://api.data.igvf.org/files/IGVFFI0000AAAA/@@upload"]


@pytest.mark.parametrize(
    "response",
    [
        _http_response(403, {"detail": "Unable to issue new credentials when upload_status..."}),
        _http_response(502, text=_GATEWAY_HTML),
    ],
)
def test_regenerate_aws_upload_creds_errors(
    monkeypatch: pytest.MonkeyPatch, connection: PConnection, response: requests.Response
) -> None:
    """Refusals and gateway errors raise HTTPError, which upload_file and its retry handle."""
    monkeypatch.setattr(connection.session, "post", lambda *_args, **_kwargs: response)
    with pytest.raises(requests.HTTPError, match=str(response.status_code)):
        connection.regenerate_aws_upload_creds("IGVFFI0000AAAA")


_PROFILES: Final[dict[str, object]] = {
    "TabularFile": {"$id": "/profiles/tabular_file.json", "properties": {}},
    "PseudobulkSet": {"$id": "/profiles/pseudobulk_set.json", "properties": {}},
    "_subtypes": {"File": ["TabularFile"]},
    "TestingDependencies": {"$id": "/profiles/testing_dependencies.json"},
    "@type": ["JSONSchemas"],
}


def test_profiles_fetched_through_session(
    monkeypatch: pytest.MonkeyPatch, connection: PConnection
) -> None:
    urls: list[str] = []

    def fake_get(url: str, **_kwargs: object) -> requests.Response:
        urls.append(url)
        return _http_response(200, _PROFILES, method="GET")

    monkeypatch.setattr(connection.session, "get", fake_get)
    assert isinstance(connection.profiles, connection_module.PProfiles)
    assert set(connection.profiles.profiles) == {"tabular_file", "pseudobulk_set"}
    schema = connection.profiles.get_profile_from_id("/tabular-files/IGVFFI0000AAAA/")
    assert schema.name == "tabular_file"
    # fetched once, then cached
    _ = connection.profiles.profiles
    assert urls == ["https://api.data.igvf.org/profiles/?format=json"]


def test_profiles_fetch_error(monkeypatch: pytest.MonkeyPatch, connection: PConnection) -> None:
    monkeypatch.setattr(
        connection.session,
        "get",
        lambda *_args, **_kwargs: _http_response(502, text=_GATEWAY_HTML, method="GET"),
    )
    with pytest.raises(requests.HTTPError, match="502"):
        _ = connection.profiles.profiles


def test_register_workers_share_profiles(connection: PConnection) -> None:
    """Each register worker uses the main thread's profiles, rather than fetching its own."""
    config = RegisterConfig(igvf_mode=IgvfMode.prod, profile_id="tabular_file")
    register._init_worker(config, connection.profiles)
    assert register._THREAD_LOCAL_DATA.connection.profiles is connection.profiles
