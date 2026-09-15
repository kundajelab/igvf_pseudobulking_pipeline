import gzip
import json
import re
import zlib
from collections.abc import Generator, Iterable, Iterator, Mapping, Sequence
from contextlib import contextmanager, nullcontext
from io import StringIO
from multiprocessing.synchronize import Lock as ProcessLock
from pathlib import Path
from threading import Lock as ThreadLock
from typing import TYPE_CHECKING, Final, Literal, cast, overload

import boto3
import igvf_utils as iu
import igvf_utils.gc_storage
import igvf_utils.utils as iuu
import requests
from igvf_utils.connection import Connection
from igvf_utils.connection import IgvfMode as IuIgvfMode
from igvf_utils.exceptions import (
    AwardPropertyMissing,
    LabPropertyMissing,
    MissingAlias,
)
from igvf_utils.profiles import IgvfSchema, Profiles
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from igvf_portal import utils
from igvf_portal.constants import VERSION
from igvf_portal.enums import (
    AnalysisStep,
    IgvfMode,
)
from igvf_portal.parallel_logger import ParallelLogger
from igvf_portal.types import (
    AccessionId,
    Alias,
    IgvfRecord,
    PortalId,
    RecordNotFound,
)

if TYPE_CHECKING:
    # type stubs only: boto3-stubs is a dev dependency, so it may be absent at runtime
    from mypy_boto3_s3 import S3Client

    # requests only defines its JSON type for type checkers
    from requests._types import JsonType

# An IGVF accession ("IGVF", two letters for the record type, four digits, four letters) names
# exactly one record, and an accessioned record's @id is /<collection>/<accession>/. So the
# accession can be read straight off either form, with no request to the Portal.
_ACCESSION_PATTERN = r"IGVF[A-Z]{2}\d{4}[A-Z]{4}"

# The Portal's web server rejects a request line longer than about 4 KB with "Request Line is too
# large" (HTTP 400). Searches went through at 3,226 characters and failed at 4,251, so stay below.
_MAX_SEARCH_URL_LENGTH: Final[int] = 3500

DEFAULT_UPLOAD_RETRY: Final[utils.RetryPolicy] = utils.RetryPolicy()
"""Default policy for retrying a failed file upload."""
_ACCESSION = re.compile(_ACCESSION_PATTERN)
_ACCESSIONED_PATH = re.compile(rf"(?:/[a-z-]+)?/({_ACCESSION_PATTERN})/")


def _local_accession(value: str) -> str | None:
    """Return the accession a value names, if it can be read off without asking the Portal.

    Handles a bare accession and a path ending in one, such as /tissues/IGVFSM0104SYAP/ or
    /IGVFSM0104SYAP/. Aliases, UUIDs, and paths of records without accessions return None.
    """
    if _ACCESSION.fullmatch(value):
        return value
    match = _ACCESSIONED_PATH.fullmatch(value)
    return match.group(1) if match else None


def _query_values(value: object) -> list[str]:
    """Spell a query parameter's value(s) as the Portal expects: a list gives one entry per item."""
    if value is None:
        values: list[object] = []
    elif isinstance(value, Sequence) and not isinstance(value, str):
        values = list(value)
    else:
        values = [value]
    # the Portal spells booleans in lower case, where str() would give "True"
    return [str(v).lower() if isinstance(v, bool) else str(v) for v in values]


_MAX_DETAIL_LENGTH: Final[int] = 1000
"""Most characters of a non-JSON response body to log."""


def _response_detail(response: requests.Response) -> str:
    """Describe a response's body for logging: formatted JSON, or else the start of its text.

    Error responses are not always JSON: a gateway in front of the Portal answers with HTML.
    """
    try:
        return iuu.print_format_dict(response.json())
    except requests.exceptions.JSONDecodeError:
        return f"HTTP {response.status_code}: {response.text[:_MAX_DETAIL_LENGTH]}"


class PProfiles(Profiles):
    """igvf_utils' Profiles, fetched through a session so that the fetch is retried."""

    def __init__(self, igvf_url: str, session: requests.Session) -> None:
        """Create profiles to be fetched from the Portal at igvf_url when first needed.

        Args:
            igvf_url: Base URL of the Portal.
            session: Session to fetch the profiles with.
        """
        super().__init__(igvf_url)
        self.session = session

    def _get_profiles(self) -> dict[str, IgvfSchema]:
        """Fetch all the public profiles on the Portal, keyed by profile ID.

        This repeats igvf_utils' Profiles._get_profiles, except for the request itself.

        Returns:
            Map from profile ID (e.g. "tabular_file") to its schema.
        """
        url = iuu.url_join([self.igvf_url, iu.PROFILES_URL, "?format=json"])
        response = self.session.get(url, timeout=iu.TIMEOUT, headers=iuu.REQUEST_HEADERS_JSON)
        response.raise_for_status()
        profiles = cast(dict[str, dict[str, object]], response.json())
        # Remove the "private" profiles (i.e. _subtypes), which have differing semantics, and the
        # "@type" pseudo profile.
        return {
            (profile_id := str(schema["$id"]).split("/")[-1].split(".json")[0]): IgvfSchema(
                profile_id, schema
            )
            for name, schema in profiles.items()
            if not (name.startswith(("_", "Testing")) or name == "@type")
        }


def _upload_error_is_transient(exception: Exception) -> bool:
    """Get whether a failed upload could succeed if tried again.

    A missing local file, or an HTTP error the Portal answers with 4xx (other than 429 Too Many
    Requests), will fail the same way every time. A 403 for a finalized file never gets here:
    upload_file skips that file.
    """
    if isinstance(exception, (FileNotFoundError, IsADirectoryError)):
        return False
    if isinstance(exception, requests.HTTPError) and exception.response is not None:
        status_code = exception.response.status_code
        return status_code >= 500 or status_code == requests.codes.TOO_MANY_REQUESTS
    return True


def _gunzip_chunks(chunks: Iterable[bytes]) -> Iterator[bytes]:
    """Incrementally gunzip a stream of chunks.

    A chunk boundary can fall anywhere in the compressed data, so chunks cannot be decompressed one
    at a time. The data can also be several gzip members end to end (as bgzip writes), so a new
    decompressor starts whenever one member ends.
    """
    decompressor = zlib.decompressobj(wbits=zlib.MAX_WBITS | 16)
    for chunk in chunks:
        while len(chunk) > 0:
            if decompressed := decompressor.decompress(chunk):
                yield decompressed
            if decompressor.eof:
                chunk = decompressor.unused_data
                decompressor = zlib.decompressobj(wbits=zlib.MAX_WBITS | 16)
            else:
                chunk = b""
    if tail := decompressor.flush():
        yield tail


class PConnection(Connection):
    """A Portal connection with retrying requests, a record cache, and S3 uploads.

    Extends igvf_utils' Connection: GET, POST, and searches go through a session that retries
    common network problems, looked-up records are cached, and a POST that conflicts with an
    existing record becomes a PUT of the differences.
    """

    session: requests.Session
    logger: ParallelLogger
    continue_on_failed_credentials: bool
    region_name: str
    record_lookups: dict[
        tuple[AccessionId | Alias | PortalId, str | None, bool], IgvfRecord | tuple[()]
    ]
    _ids_for_compare: dict[str, str]
    mode: IgvfMode

    def __init__(
        self,
        igvf_mode: IgvfMode | str,
        submission: bool = False,
        dry_run: bool = False,
        lock: ThreadLock | ProcessLock | nullcontext | None = None,
        continue_on_failed_credentials: bool = True,
        region_name: str = "us-west-2",
        upload_retry: utils.RetryPolicy = DEFAULT_UPLOAD_RETRY,
    ) -> None:
        """Connect to the Portal.

        Args:
            igvf_mode: The Portal instance to connect to, or the name of its IgvfMode.
            submission: If True, then submission mode is on, so records can be changed and GETs
                read from the database rather than the search index.
            dry_run: If True, then don't change anything on the Portal or upload any files.
            lock: Lock shared by parallel workers, used to keep log messages from interleaving.
            continue_on_failed_credentials: If True, then skip uploading a file whose upload
                credentials can't be obtained (it is probably finalized), rather than raising.
            region_name: AWS region of the S3 bucket that files are uploaded to.
            upload_retry: How to retry a failed file upload, including getting its credentials.
        """
        _igvf_mode = igvf_mode if isinstance(igvf_mode, IgvfMode) else IgvfMode[igvf_mode]
        # Pass the URL rather than the mode name: igvf_utils maps the names to its own URLs, which
        # need not agree with IgvfMode.url. Then methods inherited from Connection (profiles, patch,
        # upload credentials) could use a different Portal from the ones defined here.
        super().__init__(
            igvf_mode=_igvf_mode.url,
            submission=submission,
            dry_run=dry_run,
            no_log_file=True,
        )

        # create a session that retries common network problems
        retry = Retry(
            total=3,
            read=3,  # retries on read timeout
            connect=3,  # retries on connection timeout
            status=3,
            backoff_factor=5,  # urllib3 2 waits 0 s, 10 s, then 20 s before successive retries
            status_forcelist=[429, 500, 502, 503, 504],
            allowed_methods=["GET", "POST", "PUT", "HEAD"],  # urllib3 ≥ 1.26
            raise_on_status=False,
        )
        session = requests.Session()
        session.mount("https://", HTTPAdapter(max_retries=retry))
        session.mount("http://", HTTPAdapter(max_retries=retry))

        self.region_name = region_name
        self.session = session
        self.session.auth = self.auth
        self.logger = ParallelLogger.new(logger=self.debug_logger, lock=lock)
        self.continue_on_failed_credentials = continue_on_failed_credentials
        self.upload_retry = upload_retry
        self.record_lookups = {}
        self._ids_for_compare = {}
        self.mode = _igvf_mode

    @classmethod
    def new(
        cls,
        igvf_mode: IgvfMode,
        submission: bool = False,
        dry_run: bool = False,
        lock: ThreadLock | ProcessLock | nullcontext | None = None,
        continue_on_failed_credentials: bool = True,
        region_name: str = "us-west-2",
        upload_retry: utils.RetryPolicy = DEFAULT_UPLOAD_RETRY,
    ) -> PConnection:
        """Create a connection, fixing up logging if it runs in a process-parallel environment.

        Args:
            igvf_mode: The Portal instance to connect to.
            submission: If True, then submission mode is on, so records can be changed and GETs
                read from the database rather than the search index.
            dry_run: If True, then don't change anything on the Portal or upload any files.
            lock: Lock shared by parallel workers. A process lock means the root logger may need
                to be cleaned up again.
            continue_on_failed_credentials: If True, then skip uploading a file whose upload
                credentials can't be obtained (it is probably finalized), rather than raising.
            region_name: AWS region of the Portal's S3 backing, for public download URLs and
                uploads.
            upload_retry: How to retry a failed file upload, including getting its credentials.

        Returns:
            The new connection.
        """
        connection = cls(
            igvf_mode=igvf_mode,
            submission=submission,
            dry_run=dry_run,
            continue_on_failed_credentials=continue_on_failed_credentials,
            region_name=region_name,
            upload_retry=upload_retry,
            lock=lock,
        )
        if isinstance(lock, ProcessLock):
            # in a process-parallel environment, may need to clean up the
            # root logger again
            utils.fix_igvf_logging()

        return connection

    @property
    def igvf_mode(self) -> IuIgvfMode:
        """The Portal for methods inherited from igvf_utils.

        Overrides igvf_utils, which first checks that the Portal answers. That is unnecessary: the
        URL comes from an IgvfMode, so it is always one of the Portal's own.
        """
        return IuIgvfMode(self.mode.url)

    @property
    def profiles(self) -> Profiles:
        """The Portal's profiles, fetched through the session when first needed."""
        if self._profiles is None:
            self._profiles = PProfiles(self.mode.url, session=self.session)
        return self._profiles

    def use_profiles(self, profiles: Profiles) -> None:
        """Use profiles already fetched by another connection, rather than fetching them again.

        Profiles are only read once fetched, so connections in parallel workers can share them.
        """
        self._profiles = profiles

    def regenerate_aws_upload_creds(self, file_id: str) -> dict[str, str]:
        """Get new AWS S3 upload credentials for a file record.

        Overrides igvf_utils to go through the session, so the request is retried.

        Args:
            file_id: An identifier for a file record on the Portal.

        Returns:
            The upload_credentials of the file record.

        Raises:
            requests.exceptions.HTTPError: The Portal refused (403 if the file is finalized).
        """
        self.logger.debug(f"Getting new upload credentials for '{file_id}'")
        response = self.session.post(
            iuu.url_join([self.mode.url, "files", file_id, "@@upload"]),
            headers=iuu.REQUEST_HEADERS_JSON,
            json={},
            timeout=iu.TIMEOUT,
        )
        if not response.ok:
            self.logger.debug(
                f"Unable to get upload credentials for '{file_id}':\n{_response_detail(response)}"
            )
            response.raise_for_status()
        return cast(dict[str, str], response.json()["@graph"][0]["upload_credentials"])

    def set_submission(self, status: bool) -> None:
        """Turn submission mode on or off.

        Args:
            status: If True, then submission mode is on.
        """
        self.submission = status

    def check_dry_run(self) -> bool:
        """Return whether the dry-run feature is enabled."""
        return self.dry_run

    def _id_for_compare(self, value: str) -> str:
        """Reduce an identifier to a key that is equal exactly when two identifiers name one record.

        The key is the record's accession if it has one, otherwise its @id. An accession or a path
        ending in one is reduced locally; aliases, UUIDs, and paths of records without accessions
        are looked up once and the answer cached. Answers go in their own cache rather than
        record_lookups, because lookup_record callers expect the embedded fields that frame=object
        omits.
        """
        accession = _local_accession(value)
        if accession is not None:
            return accession
        cached = self._ids_for_compare.get(value)
        if cached is not None:
            return cached
        # Ask the search index first: it is fast, and safe here because a record's @id never
        # changes.
        # Building a record from the database can take over a minute for a heavily linked one (a
        # tissue linked by many file sets), so only use it for a record the index does not have
        # yet, such as one created moments ago.
        try:
            # try to look up via the index (much faster)
            record = self.lookup_record(cast(Alias, value), frame="object", database=False)
        except RecordNotFound:
            # try to look up via the database (if the object has just posted, it might only be
            # there)
            record = self.lookup_record(cast(Alias, value), frame="object", database=True)

        record_id = record["@id"]
        record_accession = record["accession"]
        self._ids_for_compare[value] = record_accession
        self._ids_for_compare[record_id] = record_accession
        return record_accession

    def _named_individual_props_equal[T](self, key: str, prop1: T, prop2: T) -> bool:
        if prop1 == prop2:
            return True
        if key not in {
            "file_set",
            "analysis_step_version",
            "derived_from",
            "file_format_specifications",
            "samples",
            "input_file_sets",
        }:
            return False
        if not (isinstance(prop1, str) and isinstance(prop2, str)):
            return False
        # maybe we need to look these up and compare IDs, as opposed to comparing an alias to an
        # accession
        try:
            return self._id_for_compare(prop1) == self._id_for_compare(prop2)
        except RecordNotFound:
            return False

    def _lookup_id_for_compare[T](self, key: str, value: T) -> T | str:
        lookup_props = frozenset(
            {
                "file_set",
                "analysis_step_version",
                "derived_from",
                "file_format_specifications",
                "samples",
                "input_file_sets",
            }
        )
        match value:
            case dict(d_value):
                value_id = d_value.get("@id", None)
                if value_id is not None:
                    # an embedded record's @id is canonical, so reduce it locally without a lookup
                    value_id = cast(str, value_id)
                    return _local_accession(value_id) or value_id
                return value
            case str(s_value):
                if key in lookup_props:
                    # Only a record that genuinely does not exist falls back to the raw value. Any
                    # other failure (timeout, server error) propagates: treating it as "different"
                    # would PATCH a record that did not need it.
                    try:
                        return self._id_for_compare(s_value)
                    except RecordNotFound:
                        pass
                return value
            case _:
                return value

    def _props_equal[T](self, key: str, prop1: T, prop2: T) -> bool:
        if prop1 == prop2:
            return True
        p1 = (
            sorted(self._lookup_id_for_compare(key, _s) for _s in prop1.split(","))
            if isinstance(prop1, str)
            else sorted(self._lookup_id_for_compare(key, _x) for _x in prop1)
            if isinstance(prop1, list)
            else self._lookup_id_for_compare(key, prop1)
        )
        p2 = (
            sorted(self._lookup_id_for_compare(key, _s) for _s in prop2.split(","))
            if isinstance(prop2, str)
            else sorted(self._lookup_id_for_compare(key, _x) for _x in prop2)
            if isinstance(prop2, list)
            else self._lookup_id_for_compare(key, prop2)
        )
        return type(p1) is type(p2) and p1 == p2

    def _patch_in_post(
        self,
        record_id: str,
        profile: IgvfSchema,
        payload: dict[str, object],
        existing_record: dict[str, object],
    ) -> tuple[dict[str, object], int | None]:
        """Change attempted POST that had a conflict into a patch."""
        original_record, changed, changed_props, added_props, removed_props = (
            self._get_record_changes(
                payload=payload, existing_record=existing_record, profile=profile
            )
        )

        if not changed:
            self.logger.info(f"No changes to existing record for '{record_id}', leaving as-is.")
            return original_record, None

        url = iuu.url_join([self.mode.url, record_id.lstrip("/")])
        response = self.session.put(
            url,
            timeout=iu.TIMEOUT,
            headers=iuu.REQUEST_HEADERS_JSON,
            json=cast("JsonType", existing_record),
            verify=True,
        )

        lines: list[str] = []
        if len(removed_props) > 0:
            lines.append(f"\tremoved: {','.join(prop for prop in removed_props)}")
        if len(added_props) > 0:
            lines.append(
                f"\tadded: {','.join(f'{prop}={existing_record[prop]}' for prop in added_props)}"
            )
        if len(changed_props) > 0:
            changes = ",".join(
                f"{prop}:{old_val}->{existing_record[prop]}"
                for prop, old_val in changed_props.items()
            )
            lines.append(f"changed: {changes}")
        if response.ok:
            lines.insert(0, f"Successfully PUT {record_id}.")
            self.logger.info("\n".join(lines))
        else:
            lines.insert(0, f"Failed to PUT {record_id}.")
            lines.append(_response_detail(response))
            self.logger.error("\n".join(lines))
            response.raise_for_status()
        return original_record, response.status_code

    def _get_record_changes(
        self, payload: dict[str, object], existing_record: dict[str, object], profile: IgvfSchema
    ) -> tuple[dict[str, object], bool, dict[str, object], set[str], set[str]]:
        """Find the changes in the existing record."""
        changed = False
        changed_props = {}
        added_props = set()
        removed_props = set()
        original_record = {**existing_record}
        for key, value in payload.items():
            if key in {self.IGVFID_KEY, self.PROFILE_KEY}:
                continue
            if key in ("md5sum", "file_size") and existing_record.get("upload_status") in (
                "validated",
                "validation exempted",
            ):
                # these can't be changed because the file can't be re-uploaded
                continue
            prop = profile.get_property_from_name(key)
            if existing_record.get(key) is None:
                if not prop.is_not_submittable:
                    changed = True
                    added_props.add(key)
                    existing_record[key] = value
            elif not self._props_equal(key, existing_record[key], value) and not (
                prop.is_not_submittable or prop.is_read_only
            ):
                changed = True
                changed_props[key] = existing_record[key]
                existing_record[key] = value
        to_remove_props = set(existing_record.keys()).difference(payload.keys())
        for key in to_remove_props:
            prop = profile.get_property_from_name(key)
            if prop.is_required or prop.is_not_submittable or prop.is_read_only:
                continue  # we don't remove this property
            # Then it is safe to remove this property.
            changed = True
            existing_record.pop(key)
            removed_props.add(key)
        return original_record, changed, changed_props, added_props, removed_props

    def post(
        self,
        payload: dict[str, object],
        require_aliases: bool = True,
        upload_file: bool = True,
        return_original_status_code: bool = False,
        truncate_long_strings_in_payload_log: bool = False,
        upload_duplicate: bool = True,
        expect_patch: bool = False,
    ) -> tuple[dict[str, object], int | None] | dict[str, object] | None:
        """POST a record to the Portal.

        Requires that you include in the payload the non-schematic key ``self.PROFILE_KEY`` to
        designate the name of the IGVF object profile that you are submitting to, or the
        actual `@id` property itself.

        If the `lab` property isn't present in the payload, then the default will be set to the
        value of the `IGVF_LAB` environment variable. Similarly, if the `award` property isn't
        present, then the default will be set to the value of the `IGVF_AWARD` environment
        variable.

        Before the POST is attempted, any pre-POST hooks are fist called; see the method
        ``self.before_submit_hooks``).  After a successfuly POST, any after-POST submit hooks are
        also run; see the method ``self.after_submit_hooks``.

        Args:
            payload: `dict`. The data to submit.
            require_aliases: `bool`. `True` means that the 'aliases' property is to be required in
                `payload`. This is the default and it is highly recommended not to change this
                because it'll be easy to create duplicates on the server if accidentally POSTING
                the same payload again. For example, you can easily create the same biosample
                as many times as you want on the Portal when not providing an alias. Furthermore,
                submitting labs should include at least one alias per record being submitted
                to the Portal for traceabilty purposes in the submitting lab.
            upload_file: `bool`. If `False`, when POSTing files the file data will not
                be uploaded to S3, defaults to `True`. This can be useful if you have
                custom upload logic. If the files to upload are already on disk, it is
                recommmended to leave this with the default, which will use `aws s3 cp`
                to upload them.
            return_original_status_code: `bool`. Defaults to `False`. If `True`, then
                will return the original `requests.Response.status_code` of the initial
                post, in addition to the usual `dict` response.
            truncate_long_strings_in_payload_log: `bool`. Defaults to `False`. If
                `True`, then long strings (> 1000 characters) present in the payload
                will be truncated before being logged.
            upload_duplicate: If True, upload file when a duplicate record exists. Used when
                errors allowed POSTing the record but not uploading the file.
            expect_patch: If True, then check for an existing record before attempting to post, and
                skip immediately to creating a patch if it does.

        Returns:
            `dict`: The JSON response from the POST operation, or the existing record if it already
            exists on the Portal (where a GET on any of it's aliases, when provided in the payload,
            finds the existing record). If `return_original_status_code=True`, then will
            return a `tuple` of the above `dict` and an `int` corresponding to the
            status code on POST of the initial payload.

        Raises:
            igvf_utils.exceptions.AwardPropertyMissing: The `award` property isn't present in the
                payload and there isn't a default set by the environment variable `IGVF_AWARD`.
            igvf_utils.exceptions.LabPropertyMissing: The `lab` property isn't present in the
                payload and there isn't a default set by the environment variable `IGVF_LAB`.
            igvf_utils.exceptions.MissingAlias: The argument 'require_aliases' is set to True and
                the 'aliases' property is missing in the payload or is empty.
            requests.exceptions.HTTPError: The return status is not ok.

        Side effects:
            self.PROFILE_KEY will be popped out of the payload if present, otherwise, the key "@id"
            will be popped out. Furthermore, self.IGVFID_KEY will be popped out if present in the
            payload.
        """
        self.logger.debug("\nIN post().")

        payload, profile, no_alias, aliases = self._prep_and_validate_payload(
            payload=payload, require_aliases=require_aliases
        )
        url = iuu.url_join([self.mode.url, profile.name])

        self.logger.debug(
            f"POST {profile.name} record {aliases[0]} To IGVF database with URL {url} and this "
            "payload:\n"
            + iuu.print_format_dict(
                payload, truncate_long_strings=truncate_long_strings_in_payload_log
            )
        )

        if self.check_dry_run():
            return {}
        if expect_patch and not no_alias:
            existing_record = self.get(rec_ids=aliases, ignore404=True, frame="edit")
            if existing_record is not None:
                # try to immediately go for a patch
                return self._handle_conflict(
                    payload=payload,
                    existing_record=existing_record,
                    aliases=aliases,
                    profile=profile,
                    upload_file=upload_file,
                    upload_duplicate=upload_duplicate,
                    return_original_status_code=return_original_status_code,
                )

        # try to post de-novo
        response = self.session.post(
            url,
            timeout=iu.TIMEOUT,
            headers=iuu.REQUEST_HEADERS_JSON,
            json=cast("JsonType", payload),
            verify=True,
        )
        original_status_code = response.status_code

        # NOTE: only parse the body as JSON once the POST is known to have succeeded. An error can
        # come from a gateway as HTML, and parsing it would hide the status code.
        if response.ok:
            self.logger.debug("Success.")
            response_json = response.json()["@graph"][0]
            # Some objects don't have an accession, i.e. replicates.
            # in original code "record_id" was frequently called "encid" (presumably encode id?)
            record_id: str = AccessionId(response_json.get("accession", response_json["uuid"]))
            self.logger.debug(f"Object posted with identifier: {record_id}")
            self._log_post(aliases=aliases, dacc_id=record_id)
            # Run 'after' hooks:
            self.guard_upload(
                upload_file=upload_file,
                profile=profile,
                accession_id=record_id,
                payload=payload,
            )
            if return_original_status_code is True:
                return (response_json, original_status_code)
            return response_json
        elif response.status_code == requests.codes.CONFLICT:
            # In the case of paired-end FASTQ files, it could also mean that there was a conflict
            # related to the 'paired_with' property, i.e. the latter is already linked to a FASTQ
            # file, which could even have been set to a deleted state on the Portal. The server
            # response in either case would look something like this:
            #
            # {
            #   'detail':
            #     "Keys conflict: [('file:paired_with', 'f39320d9-0970-4369-b680-5965a5e85b6f')]",
            #   'description': 'There was a conflict when trying to complete your request.',
            #   'code': 409,
            #   '@type': ['HTTPConflict', 'Error'],
            #   'title': 'Conflict',
            #   'status': 'error'}
            # }
            #
            if no_alias:
                self.logger.warning(_response_detail(response))
                response.raise_for_status()
            else:
                existing_record = self.get(rec_ids=aliases, ignore404=True, frame="edit")
                if existing_record is None:
                    self.logger.warning(_response_detail(response))
                    response.raise_for_status()
                else:
                    return self._handle_conflict(
                        payload=payload,
                        existing_record=existing_record,
                        aliases=aliases,
                        profile=profile,
                        upload_file=upload_file,
                        upload_duplicate=upload_duplicate,
                        return_original_status_code=return_original_status_code,
                    )

        else:
            self.logger.error(f"Failed to POST {aliases[0]}\n{_response_detail(response)}")
            response.raise_for_status()

    def _prep_and_validate_payload(
        self,
        payload: dict[str, object],
        require_aliases: bool,
    ) -> tuple[dict[str, object], IgvfSchema, bool, list[Alias]]:
        # Make sure we have a payload that can be converted to valid JSON, and
        # tuples become arrays, ...
        payload = json.loads(json.dumps(payload))
        profile = self.get_profile_from_payload(payload)
        payload[self.PROFILE_KEY] = profile.name
        # Check if we need to add defaults for 'award' and 'lab' properties:
        if profile.has_award:  # No lab prop for these profiles either.
            if iu.AWARD_PROP_NAME not in payload:
                if not iu.AWARD:
                    raise AwardPropertyMissing
                payload.update(iu.AWARD)
            if iu.LAB_PROP_NAME not in payload:
                if not iu.LAB:
                    raise LabPropertyMissing
                payload.update(iu.LAB)

        # Run 'before' hooks:
        payload = self.before_submit_hooks(payload, method=self.POST)

        # Remove the non-schematic self.PROFILE_KEY if being used, which was added above since some
        # 'before' hooks may need it. Also check for the `@id` property and remove it too if found.
        def _is_wanted_key(_k: str) -> bool:
            if _k in {self.PROFILE_KEY, "@id", self.IGVFID_KEY}:
                return False
            _prop = profile.get_property_from_name(_k)
            return not (_prop.is_not_submittable or _prop.is_read_only)

        payload = {k: v for k, v in payload.items() if _is_wanted_key(k)}

        no_alias, aliases = self._get_aliases(
            payload=payload, profile=profile, require_aliases=require_aliases
        )

        # Validate the payload against the schema
        ### This doesn't work as locally I can't use jsonschema to validate a profile with
        ### custom objects specified in the value of a linkTo property.
        self.logger.debug("Validating the payload against the schema")
        validation_error = iuu.err_context(
            payload=payload,
            schema=self.profiles.get_profile_from_id(profile.name).schema,
        )
        if validation_error:
            self.logger.error(f"Invalid schema instance of the {profile.name} profile.")
            self.logger.error(f"Payload is: {iuu.print_format_dict(payload)}")
            self.logger.error(validation_error[0])  # The top-level validation message
            if validation_error[1]:  # The validation context can be empty
                self.logger.error(iuu.print_format_dict(validation_error[1]))
            raise Exception(iuu.print_format_dict(validation_error[0]))

        return payload, profile, no_alias, aliases

    @classmethod
    def _get_aliases(
        cls, payload: dict[str, object], profile: IgvfSchema, require_aliases: bool = True
    ) -> tuple[bool, list[Alias]]:
        """Get the list of Aliases, and a boolean if the original payload is missing one."""
        no_alias = False  # Use this to check later if doing a GET
        aliases = cast(list[Alias], payload.get(iu.ALIAS_PROP_NAME))
        if not aliases:
            if not profile.has_alias or not require_aliases:
                aliases = [Alias("N/A")]
                no_alias = True
            else:
                raise MissingAlias(
                    f"Missing property '{iu.ALIAS_PROP_NAME}' in payload {payload}. This is"
                    " required by default for the profiles that include this property, and can be"
                    " disabled by setting the `require_aliases` argument to False in the call to"
                    " this method, being `igvf_utils.connection.Connection.post()`. If using the"
                    " iu_register.py script, this can also be disabled by passing the"
                    " --no-aliases option."
                )
        return no_alias, aliases

    def _handle_conflict(
        self,
        payload: dict[str, object],
        existing_record: dict[str, object],
        aliases: list[Alias],
        profile: IgvfSchema,
        upload_file: bool = True,
        upload_duplicate: bool = True,
        return_original_status_code: bool = False,
    ) -> tuple[dict[str, object], int | None] | dict[str, object]:
        if upload_duplicate:
            accession_id = cast(AccessionId, existing_record.get("accession", aliases[0]))
            record_id = cast(str, existing_record.get("@id", accession_id))
            self.logger.warning(f"Conflict when POSTing {record_id}, will patch any differences.")
            original_record, status_code = self._patch_in_post(
                record_id=record_id,
                profile=profile,
                payload=payload,
                existing_record=existing_record,
            )
            self.guard_upload(
                upload_file=upload_file,
                profile=profile,
                accession_id=accession_id,
                payload=payload,
                original_record=original_record,
            )
            return (
                (original_record, status_code) if return_original_status_code else original_record
            )
        else:
            self.logger.error(
                f"Will not POST '{aliases[0]}' since it already exists with aliases "
                f"'{existing_record['aliases']}'."
            )
        return (existing_record, None) if return_original_status_code else existing_record

    def _regenerate_s3_client(self, file_id: str) -> tuple[str, str, S3Client]:
        """Get a client for uploading the specified client to S3.

        Args:
            file_id: File identifier (PortalId, Alias, or Accession)

        Returns:
            bucket: str with S3 bucket name
            key: str with path to file_id in bucket
            s3_client: client to interact with S3
        Throws:
            requests.exceptions.HTTPError if credentials cannot be obtained
        """
        upload_credentials = self.regenerate_aws_upload_creds(file_id)
        aws_creds = {
            "AWS_ACCESS_KEY_ID": upload_credentials["access_key"],
            "AWS_SECRET_ACCESS_KEY": upload_credentials["secret_key"],
            "AWS_SESSION_TOKEN": upload_credentials["session_token"],
            "UPLOAD_URL": upload_credentials["upload_url"],
            "AWS_ACCOUNT_ID": upload_credentials["federated_user_arn"],
        }
        s3_client = boto3.client(
            "s3",
            region_name=self.region_name,
            aws_access_key_id=aws_creds["AWS_ACCESS_KEY_ID"],
            aws_secret_access_key=aws_creds["AWS_SECRET_ACCESS_KEY"],
            aws_session_token=aws_creds["AWS_SESSION_TOKEN"],
            aws_account_id=aws_creds["AWS_ACCOUNT_ID"],
        )
        bucket, key = utils.parse_s3_uri(aws_creds["UPLOAD_URL"])
        return bucket, key, s3_client

    def _is_already_uploaded(
        self,
        file_id: str,
        bucket: str,
        key: str,
        s3_client: S3Client,
        local_md5_sum: str,
    ) -> bool:
        # subprocess_result = subprocess.run(
        #     f"aws s3api head-object --bucket '{bucket}' --key '{key}'",
        #     shell=True,
        #     capture_output=True,
        #     env={**os.environ, **aws_creds},
        # )
        # if subprocess_result.returncode != 0:
        #     self.logger.info(f"Upload url '{upload_url}' doesn't already exist.")
        #     return False
        # result = json.loads(subprocess_result.stdout)

        self.logger.debug(f"Getting head info for Bucket={bucket}, Key={key}")
        try:
            result = s3_client.head_object(Bucket=bucket, Key=key)
        except s3_client.exceptions.NoSuchKey, s3_client.exceptions.ClientError:
            self.logger.debug(
                f"Upload url 's3://{bucket}/{key}' doesn't already exist for '{file_id}'."
            )
            return False
        except Exception as exception:
            raise RuntimeError(
                f"Getting head info for {file_id} at 's3://{bucket}/{key}'"
            ) from exception

        remote_md5_sum = result.get("Metadata", {}).get("md5sum", "")
        # this is unneccessary: interrupted/failed uploads do not result in remote obects,
        # so the only way the metadata could be present is if the object uploads successfully
        # remote_crc64_nvme_checksum = result.get("ChecksumCRC64NVME", "")
        # local_crc64_nvme_checksum = utils.crc64_nvme_checksum(file_path)
        already_uploaded = remote_md5_sum == local_md5_sum
        if already_uploaded:
            self.logger.info(
                f"Upload url 's3://{bucket}/{key}' md5_sum matches local file, will not re-upload."
            )
        else:
            self.logger.info(
                f"Upload url 's3://{bucket}/{key}' md5_sum does NOT match local file, will "
                "re-upload."
            )
        return already_uploaded

    def guard_upload(
        self,
        upload_file: bool,
        profile: IgvfSchema,
        accession_id: AccessionId,
        payload: dict[str, object],
        original_record: dict[str, object] | None = None,
    ) -> None:
        """Upload the file for a file record, if it should be uploaded.

        Nothing is uploaded if upload_file is False, the record is not a file, or the existing
        record's file is already validated (or exempted from validation).

        Args:
            upload_file: If False, then don't upload anything.
            profile: Profile of the record.
            accession_id: Accession of the record.
            payload: The record's submitted payload, which names the local file to upload.
            original_record: The record as it was on the Portal before this submission, if it
                already existed.
        """
        if upload_file:
            if profile.name in self.profiles.FILE_PROFILE_ID:
                file_path = Path(str(payload[self.profiles.SUBMITTED_FILE_PROP_NAME]))
                md5sum = payload.get("md5sum")
                set_md5sum = md5sum if isinstance(md5sum, str) and len(md5sum) > 0 else None
                if original_record is not None:
                    upload_status = cast(str | None, original_record.get("upload_status", None))
                    if set_md5sum is None:
                        set_md5sum = utils.md5sum(file_path)
                    original_md5sum = cast(str | None, original_record.get("md5sum", None))
                    if upload_status in ("validated", "validation exempted"):
                        if original_md5sum is not None and original_md5sum != set_md5sum:
                            self.logger.warning(
                                f"Skipping upload of {accession_id} with changed md5sum because "
                                f"remote file is {upload_status}. Contact the DACC to invalidate "
                                "if it needs to be replaced."
                            )
                        else:
                            self.logger.info(
                                f"Skipping upload of {accession_id} because remote file is "
                                f"{upload_status}."
                            )
                        return
                    else:
                        if original_md5sum is not None and original_md5sum != set_md5sum:
                            self.logger.info(f"Uploading {accession_id} with changed md5sum.")
                        else:
                            self.logger.info(
                                f"Uploading {accession_id} with unchanged md5sum, because "
                                f"upload_status is {upload_status} and S3 backing cannot be "
                                "directly checked."
                            )
                # Retry the whole upload: each attempt gets fresh credentials, so this also recovers
                # from credentials that expire part-way through a long upload.
                utils.retry(
                    num_tries=self.upload_retry.num_tries,
                    delay=self.upload_retry.delay,
                    backoff=self.upload_retry.backoff,
                    should_retry=_upload_error_is_transient,
                    logger=self.logger.logger,
                )(self.upload_file)(
                    file_id=accession_id, file_path=file_path, set_md5sum=set_md5sum
                )
            else:
                self.logger.debug(f"No upload of {accession_id} because record is not a file.")
        else:
            self.logger.debug(f"Skipping upload of {accession_id} because upload is False")

    def upload_file(
        self,
        file_id: AccessionId,
        file_path: str | Path | None = None,
        set_md5sum: str | bool | None = None,
    ) -> None:
        """Upload a file to the Portal for the indicated file record.

        The file to upload can be specified by setting the `file_path` parameter, or by using the
        value of the IGVF file profile's `submitted_file_name` property of the given file object
        represented by the `file_id` parameter. The file to upload can be from any of the following
        sources:

          1. Path to a local file,
          2. S3 object, or
          3. Google Storage object

        For the AWS option above, the user must set the proper AWS keys, see the
        `wiki documentation`_.

        If the dry-run feature is enabled, then this method will return prior to launching the
        upload command.

        Args:
            file_id: `str`. An identifier of a `file` record on the IGVF Portal.
            file_path: `str`. The local path to the file to upload, or an S3 object (i.e
              s3://mybucket/test.txt), or a Google Storage object (i.e. gs://mybucket/test.txt).
              If not set, defaults to `None` in which case the local file path will be extracted
              from the record's `submitted_file_name` property.
            set_md5sum: Ignored unless it is a string. Then if non-empty it is the value of the
              file's local md5sum.

        Raises:
            igvf_utils.exceptions.FileUploadFailed: The return code of the AWS upload command was
              non-zero.

        .. _`wiki documentation`: https://github.com/IGVF-DACC/igvf_utils/wiki/Configuration#aws-keys
        """
        self.logger.debug("\nIN upload_file()\n")
        # upload_credentials = self.get_upload_credentials(file_id)
        # Don't use this - they may have expired.
        try:
            bucket, key, s3_client = self._regenerate_s3_client(file_id)
        except requests.exceptions.HTTPError as http_error:
            # The Portal refuses new credentials with 403 Forbidden when the file is finalized
            # (validated or validation exempted, no longer in progress, or externally hosted).
            # Any other failure, e.g. a server error, means the file still needs uploading.
            forbidden = (
                http_error.response is not None
                and http_error.response.status_code == requests.codes.FORBIDDEN
            )
            if self.continue_on_failed_credentials and forbidden:
                self.logger.warning(
                    f"Skipping upload of {file_id} because the Portal refused upload credentials. "
                    "It is probably finalized."
                )
                return
            else:
                raise
        match file_path:
            case None:
                file_rec: IgvfRecord = self.lookup_record(file_id, frame="object")
                try:
                    file_path: Path = Path(file_rec["submitted_file_name"])
                except KeyError as key_error:  # submitted_file_name property not set:
                    raise Exception("No file path specified.") from key_error
            case str(str_path):
                file_path = Path(str_path)
            case Path():
                pass
            case _:
                raise ValueError(f"Invalid type {type(file_path)} for file_path.")

        if isinstance(set_md5sum, str) and len(set_md5sum) > 0:
            md5sum = set_md5sum
        else:
            md5sum = utils.md5sum(file_path)

        self.logger.info(f"Uploading {file_path} to 's3://{bucket}/{key}'")
        if self.check_dry_run():
            return

        s3_client.upload_file(
            Filename=f"{file_path}",
            Bucket=bucket,
            Key=key,
            ExtraArgs={"Metadata": {"md5sum": md5sum}},
        )
        self.logger.info(f"Successfully uploaded '{file_path}' to 's3://{bucket}/{key}'.")

    def get(
        self,
        rec_ids: str | Sequence[str],
        database: bool | None = None,
        ignore404: bool = True,
        frame: str | None = None,
    ) -> dict[str, object] | None:
        """GET a record from the Portal.

        Looks up a record in the Portal and performs a GET request, returning the JSON
        serialization of the object. You supply a list of identifiers for a specific record, and
        the Portal will be searched for each identifier in turn until one is either found or the
        list is exhausted.

        Args:
            rec_ids: `str` or `list`. Must be a `list` if you want to supply more than one
                identifier. For a few example identifiers, you can use a uuid, accession, ..., or
                even the value of a record's `@id` property.
            database: `bool`. If True, then search the database directly instead of the
                Elasticsearch indices. Defaults to True when in submission mode (`self.submission`
                is True), otherwise False.
            frame: `str`. A value for the frame query parameter, i.e. 'object', 'edit'.
            ignore404: `bool`. Only matters when none of the passed in record IDs were found on the
                Portal.  In this case, If set to `True`, then None will be returned.
                If set to `False`, then an Exception will be raised.


        Returns:
            `dict`: The JSON response. Will be empty if no record was found AND ``ignore404=True``.

        Raises:
            ValueError: No identifiers were given.
            `Exception`: If the server responds with a FORBIDDEN status.
            `requests.exceptions.HTTPError`: The status code is not ok, and the
                cause isn't due to a 404 (not found) status code when ``ignore404=True``.
        """
        if database is None:
            database = self.submission
        if isinstance(rec_ids, str):
            rec_ids = [rec_ids]
        status_codes = {}
        response: requests.Response | None = None
        for record_id in rec_ids:
            response = self._request_record(record_id=record_id, database=database, frame=frame)
            if response.ok:
                return response.json()
            status_codes[response.status_code] = record_id
        if response is None:
            raise ValueError("No record identifiers given.")

        if requests.codes.FORBIDDEN in status_codes:
            raise RuntimeError(
                f"Access to IGVF record {status_codes[requests.codes.FORBIDDEN]} is forbidden"
            )
        elif requests.codes.NOT_FOUND in status_codes:
            self.logger.debug("NOT FOUND")
            if ignore404:
                return None
            try:
                response.raise_for_status()
            except requests.HTTPError as http_error:
                raise RecordNotFound(f"{','.join(rec_ids)} not found") from http_error

        # At this point in the code, the response is not okay.
        # Raise the error for last response we got:
        response.raise_for_status()

    def _request_record(
        self, record_id: str, database: bool, frame: str | None
    ) -> requests.Response:
        """Request the record for an individual key."""
        record_id = record_id.strip("/")
        url = iuu.url_join([self.mode.url, record_id, "?format=json"])
        if database:
            url += "&datastore=database"
        if frame:
            url += f"&frame={frame}"
        self.logger.debug(f"GET '{record_id}' From DACC with URL '{url}'")
        return self.session.get(
            url,
            timeout=iu.TIMEOUT,
            headers=iuu.REQUEST_HEADERS_JSON,
            verify=True,
        )

    def get_record_http_url(
        self,
        file_record: IgvfRecord,
    ) -> str:
        """Get a public HTTPS URL for downloading the specified file record."""
        s3_uri = file_record.get("s3_uri", "")
        if s3_uri.startswith("s3://igvf-public/"):
            # it's a public S3 URL, construct equivalent https URL and return it
            _, _, bucket_name, object_key = s3_uri.split("/", 3)
            return f"https://{bucket_name}.s3.{self.region_name}.amazonaws.com/{object_key}"
        else:
            # it's a private S3 URL, but we can request a temporary download via the "href" field
            response = self.session.head(
                url=f"{self.mode.url}{file_record['href']}",
                allow_redirects=True,
                timeout=iu.TIMEOUT,
            )
            # The final response is not checked: the Portal redirects to an S3 URL presigned for
            # GET, so S3 refuses this HEAD with 403 even though the URL is good. Without any
            # redirect, though, the Portal itself refused (e.g. bad credentials or missing file),
            # and response.url would be the Portal's own href, which aria2c cannot download without
            # credentials.
            if len(response.history) == 0:
                response.raise_for_status()
                raise RuntimeError(
                    f"Portal did not redirect {file_record['href']} to a download URL "
                    f"(status {response.status_code})."
                )
            return response.url

    def search_records(
        self,
        record_type: str | Sequence[str] | None = None,
        *,
        query: str | None = None,
        field_filters: Mapping[str, object] | None = None,
        field: str | Sequence[str] | None = None,
        frame: str | None = "object",
        limit: int | Literal["all"] = "all",
        sort: str | Sequence[str] | None = None,
    ) -> list[dict[str, object]]:
        """Search the Portal's index and return the matching records.

        Takes the same arguments as igvf_client's IgvfApi.search (except type->record_type), and
        like it asks for frame=object by default, but returns plain dicts rather than typed models.
        A caller that wants a model can build one, e.g. AnalysisSet.from_dict(record). (Named
        search_records because the igvf_utils base class already has an unrelated search().)

        Args:
            record_type: Item type(s) to search, e.g. "AnalysisSet". Several types match any of
                them.
            query: Free-text query.
            field_filters: Filters by property, e.g. {"status!": "deleted", "files.content_type":
                "cell annotations"}. End a name with "!" to negate it; a list value matches any of
                its items. Values are passed through, so range forms such as "gte:30000" work.
            field: Return only these properties (plus @id and @type). This takes precedence over
                frame, and can make a large search far smaller.
            frame: Which view of each record to return; None leaves it to the Portal.
            limit: Maximum number of results; "all" (the default) returns every match.
            sort: Properties to sort by, prefixed with "-" for descending. Does not work with
                limit="all".

        Returns:
            The matching records, in the Portal's order. Empty if nothing matches.

        Raises:
            requests.exceptions.HTTPError: The search failed.
        """
        url = self._search_url(
            record_type=record_type,
            query=query,
            field_filters=field_filters,
            field=field,
            frame=frame,
            limit=limit,
            sort=sort,
        )
        if len(url) > _MAX_SEARCH_URL_LENGTH:
            raise ValueError(
                f"Search URL is {len(url)} characters, but the Portal rejects URLs over about"
                f" {_MAX_SEARCH_URL_LENGTH}. To match any of a long list of values, use"
                " search_records_matching_any, which splits them across several searches."
            )
        self.logger.debug(f"Search '{url}'")
        response = self.session.get(url, headers=iuu.REQUEST_HEADERS_JSON, timeout=iu.TIMEOUT)
        if response.status_code == requests.codes.NOT_FOUND:
            # The Portal answers a search that matches nothing with 404, but still sends a search
            # response with an empty @graph. A 404 without one is a real failure, e.g. a bad URL.
            try:
                body = response.json()
            except requests.exceptions.JSONDecodeError:
                body = None
            if isinstance(body, dict) and "@graph" in body:
                return []
        response.raise_for_status()
        return cast(list[dict[str, object]], response.json()["@graph"])

    def download_bytes(self, key: AccessionId) -> bytes:
        """Download the content of a file record into memory.

        GETs /<key>/@@download, which the Portal answers with a redirect to an S3 link. The
        connection's credentials are needed for unreleased files; requests drops them when the
        redirect leaves the Portal, which S3 requires.
        """
        url = f"{self.mode.url.rstrip('/')}/{str(key).strip('/')}/@@download"
        self.logger.debug(f"Download '{url}'")
        response = self.session.get(url, timeout=iu.TIMEOUT)
        response.raise_for_status()
        return response.content

    @overload
    def stream_bytes(
        self, key: AccessionId, decompress: bool, chunk_size: int
    ) -> Iterator[bytes]: ...

    @overload
    def stream_bytes(self, key: AccessionId, decompress: bool, chunk_size: None) -> bytes: ...

    def stream_bytes(
        self,
        key: AccessionId | IgvfRecord,
        decompress: bool = True,
        chunk_size: int | None = 2**20,
    ) -> bytes | Iterator[bytes]:
        """Stream the content of a file record into memory.

        GETs /<key>/@@download, which the Portal answers with a redirect to an S3 link. The
        connection's credentials are needed for unreleased files; requests drops them when the
        redirect leaves the Portal, which S3 requires.

        Args:
            key: Accession ID of record in Portal, or a record that has already been looked up.
            decompress: If True, use gzip to decompress chunks before returning.
            chunk_size: If chunk_size is None, get the entire file and return the bytes.
                If chunk_size is an integer, yield chunks of that size (or smaller).

        Returns:
            bytes of chunk_size is None, otherwise Iterator of bytes
        """
        url = (
            f"{self.mode.url.rstrip('/')}/{str(key).strip('/')}/@@download"
            if isinstance(key, str)
            else self.get_record_http_url(key)
        )

        if chunk_size is None:
            self.logger.debug(f"Download '{url}'")
            response = self.session.get(url, timeout=iu.TIMEOUT)
            response.raise_for_status()
            return gzip.decompress(response.content) if decompress else response.content
        # NOTE: the streaming is in a separate generator method: a yield here would make this whole
        # method a generator, so the chunk_size=None branch could not return bytes
        return self._iter_download_chunks(url, decompress=decompress, chunk_size=chunk_size)

    def _iter_download_chunks(self, url: str, decompress: bool, chunk_size: int) -> Iterator[bytes]:
        """Yield the content downloaded from url in chunks, gunzipping it if decompress is True."""
        self.logger.debug(f"Stream '{url}'")
        with self.session.get(url, stream=True, timeout=iu.TIMEOUT) as remote_in:
            remote_in.raise_for_status()
            chunks = remote_in.iter_content(chunk_size=chunk_size)
            yield from _gunzip_chunks(chunks) if decompress else chunks

    @contextmanager
    def read_remote_bytes(self, key: AccessionId, decompress: bool = True) -> Generator[StringIO]:
        """Open bytes from IGVF Record file as TextIO.

        Args:
            key: Accession ID of record in Portal
            decompress: If True, use gzip to decompress chunks before returning.

        Yields:
            StringIO for reading
        """
        # only wrap download failures: an exception raised in the caller's with-block must propagate
        try:
            text = self.stream_bytes(key=key, chunk_size=None, decompress=decompress).decode()
        except Exception as exception:
            raise RuntimeError(f"Could not open {key} for remote reading") from exception
        yield StringIO(text)

    def download_record(
        self,
        record: IgvfRecord,
        chunk_size: int = 2**20,
        output: Path | None = None,
    ) -> Path:
        """Download the file in the provided record, and return the output path.

        Args:
            record: IGVF portal record to download.
            chunk_size: Chunk size for streaming download, in bytes.
            output: If specified, download to that folder using href as file name.
              If unspecified, download in working folder.

        Returns:
            Final local path to downloaded file.
        """
        download_path = Path(record["href"].rsplit("/", 1)[-1])
        if output is None:
            output = download_path
        elif output.is_dir() or len(output.suffixes) == 0:
            # assume output specified the folder, not the whole path, so fix it
            output = output / download_path.name

        output.parent.mkdir(exist_ok=True, parents=True)

        with output.open("wb") as local_out:
            for chunk in self.stream_bytes(
                key=record["accession"], decompress=False, chunk_size=chunk_size
            ):
                local_out.write(chunk)
        return output

    def search_records_matching_any(
        self,
        record_type: str | Sequence[str] | None,
        filter_name: str,
        values: Iterable[str],
        *,
        field_filters: Mapping[str, object] | None = None,
        field: str | Sequence[str] | None = None,
        frame: str | None = "object",
    ) -> list[dict[str, object]]:
        """Search for records whose filter_name property matches any of values, however many.

        A repeated filter matches any of its values, but the Portal rejects a URL over about 4 KB,
        so a long list of values is split across as many searches as that needs. Each record is
        returned once, in the order first found. Other arguments are as for search_records.
        """
        unique_values = sorted(set(values))
        records: dict[object, dict[str, object]] = {}

        def filters_for(batch: list[str]) -> dict[str, object]:
            return {**(field_filters or {}), filter_name: batch}

        start = 0
        while start < len(unique_values):
            # grow the batch while its URL still fits (a value too long on its own is left to
            # search_records to report)
            end = start + 1
            while end < len(unique_values) and (
                len(
                    self._search_url(
                        record_type=record_type,
                        field_filters=filters_for(unique_values[start : end + 1]),
                        field=field,
                        frame=frame,
                    )
                )
                <= _MAX_SEARCH_URL_LENGTH
            ):
                end += 1
            for record in self.search_records(
                record_type,
                field_filters=filters_for(unique_values[start:end]),
                field=field,
                frame=frame,
            ):
                records.setdefault(record.get("@id"), record)
            start = end
        return list(records.values())

    def _search_url(
        self,
        record_type: str | Sequence[str] | None = None,
        *,
        query: str | None = None,
        field_filters: Mapping[str, object] | None = None,
        field: str | Sequence[str] | None = None,
        frame: str | None = "object",
        limit: int | Literal["all"] = "all",
        sort: str | Sequence[str] | None = None,
    ) -> str:
        """The full URL search_records sends for these arguments."""
        params: list[tuple[str, str]] = [("format", "json"), ("limit", str(limit))]
        if query is not None:
            params.append(("query", query))
        for name, value in (("type", record_type), ("sort", sort), ("field", field)):
            params.extend((name, v) for v in _query_values(value))
        for name, value in (field_filters or {}).items():
            params.extend((name, v) for v in _query_values(value))
        if frame is not None:
            params.append(("frame", frame))
        base = f"{self.mode.url.rstrip('/')}/search/"
        return cast(str, requests.Request("GET", base, params=params).prepare().url)

    def _get_record_from_lookups[T: AccessionId | Alias | PortalId](
        self, keys: Iterable[T], frame: str | None, database: bool
    ) -> IgvfRecord | tuple[()] | None:
        """Retrieve record from cached lookups, or return None if it is not present."""
        all_negative_cache = True
        for key in keys:
            match self.record_lookups.get((key, frame, database), None):
                case dict(record):
                    return record
                case tuple():
                    pass
                case None:
                    all_negative_cache = False
        return () if all_negative_cache else None

    def _cache_record[T: AccessionId | Alias | PortalId](
        self,
        keys: Iterable[T],
        frame: str | None,
        database: bool,
        record: IgvfRecord | tuple[()],
    ) -> None:
        """Store retrieved record in record_lookups."""
        keys_to_cache: Iterable[T]
        match record:
            case tuple():
                keys_to_cache = keys
            case _:
                extra_keys = [
                    record.get("accession", None),
                    record.get("@id", None),
                    *record.get("aliases", []),
                ]
                keys_to_cache = set(keys)
                keys_to_cache.update(key for key in extra_keys if key is not None)  # ty:ignore[invalid-argument-type]

        for key in keys_to_cache:
            self.record_lookups[(key, frame, database)] = record

    def lookup_record[T: AccessionId | Alias | PortalId](
        self,
        key: T | Iterable[T],
        database: bool | None = None,
        ignore404: bool = True,
        frame: str | None = None,
        cache_negative: bool = False,
    ) -> IgvfRecord:
        """Look up a record, using the table of previous lookups if it's present there.

        Otherwise get it directly from the Portal and add it to the table.

        Args:
            key: An identifier, or several identifiers, of the record. For a few example
                identifiers, you can use a uuid, accession, ..., or even the value of a record's
                `@id` property.
            database: `bool`. If True, then search the database directly instead of the
                Elasticsearch indices. Default True when in submission mode (`self.submission` is
                True), otherwise False.
            ignore404: `bool`. Only matters when none of the passed in record IDs were found on the
                Portal.  In this case, If set to `True`, then None will be returned.
                If set to `False`, then an Exception will be raised.
            frame: `str`. A value for the frame query parameter, i.e. 'object', 'edit'.
            cache_negative: If true, also store 404s and throw RecordNotFound without querying.

        Returns:
            The record.

        Raises:
            RecordNotFound: No record was found for any of the identifiers.
        """
        keys: list[T] = [cast(T, key)] if isinstance(key, str) else list(key)
        if database is None:
            database = self.submission
        match self._get_record_from_lookups(keys=keys, frame=frame, database=database):
            case dict(record):
                return record
            case ():
                raise RecordNotFound(f"Could not find record for '{','.join(keys)}'")
            case None:
                record = cast(
                    IgvfRecord | None,
                    self.get(rec_ids=keys, database=database, ignore404=ignore404, frame=frame),
                )
                if record is None:
                    if cache_negative:
                        self._cache_record(keys=keys, frame=frame, database=database, record=())
                    raise RecordNotFound(f"Could not find record for '{','.join(keys)}'")
                self._cache_record(keys=keys, frame=frame, database=database, record=record)
                return record

    def infer_principal_accessions(
        self, intermediate_accessions: Iterable[AccessionId]
    ) -> set[PortalId]:
        """Check intermediate accessions to find the principal accessions they derive from."""
        # check all the supplied intermediate accessions
        to_check = set(intermediate_accessions)
        principal_ids: set[PortalId] = set()
        while len(to_check) > 0:
            # pop off one of the intermediate accessions and get its record
            intermediate_accession = to_check.pop()
            intermediate_record = self.lookup_record(intermediate_accession)
            if intermediate_record["status"] == "deleted":
                continue
            principal_ids_for_intermediate = intermediate_record.get("input_for", None)
            if principal_ids_for_intermediate is None:
                # it's not input for anything, so it must be a principal accession
                principal_ids.add(intermediate_record["@id"])
            else:
                # it's input for these principal ids.
                for principal_id in principal_ids_for_intermediate:
                    # get the record for this principal ID
                    principal_record = self.lookup_record(principal_id)
                    if principal_record["status"] == "deleted":
                        continue
                    # get its accession and add it to the output set
                    principal_ids.add(principal_record["@id"])
                    # to decrease lookups, remove everything that was input to it from the IDs to
                    # check
                    # (this is VERY effective for data sets with many intermediate accessions)
                    intermediate_inputs = (
                        rec["accession"] for rec in principal_record.get("input_file_sets", [])
                    )
                    to_check.difference_update(intermediate_inputs)

        return principal_ids

    def lookup_analysis_step_version(self, analysis_step: AnalysisStep) -> list[Alias]:
        """Find the aliases of this pipeline version's analysis step version for an analysis step.

        Args:
            analysis_step: The analysis step.

        Returns:
            Aliases of the analysis step version whose software version is this pipeline's.

        Raises:
            ValueError: The analysis step has no version for this pipeline version.
        """
        step_record = cast(dict[str, object], self.lookup_record(analysis_step.value))
        for version_dict in cast(list[dict[str, object]], step_record["analysis_step_versions"]):
            for software_versions in cast(list[dict[str, str]], version_dict["software_versions"]):
                if software_versions["name"] == f"igvf_pseudobulking_pipeline-v{VERSION}":
                    version_id = cast(Alias, version_dict["@id"])
                    return self.lookup_record(version_id, frame="object")["aliases"]
        raise ValueError(
            f"Unable to find version of analysis step {analysis_step} for "
            f"'igvf_pseudobulking_pipeline-v{VERSION}'"
        )
