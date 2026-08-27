import json
import os
import subprocess
from contextlib import nullcontext
from multiprocessing.synchronize import Lock as ProcessLock
from pathlib import Path
from threading import Lock as ThreadLock

import igvf_utils as iu
import igvf_utils.gc_storage
import igvf_utils.transfer_to_gcp
import igvf_utils.utils as iuu
import requests
from igvf_utils.connection import Connection
from igvf_utils.exceptions import (
    AwardPropertyMissing,
    LabPropertyMissing,
    MissingAlias,
)
from igvf_utils.profiles import IgvfSchema

from igvf_portal import utils
from igvf_portal.enums import IgvfMode
from igvf_portal.parallel_logger import ParallelLogger
from igvf_portal.types import IgvfRecord


class PConnection(Connection):
    session: requests.Session
    logger: ParallelLogger
    continue_on_failed_credentials: bool

    def __init__(
        self,
        igvf_mode: IgvfMode,
        submission: bool = False,
        dry_run: bool = False,
        lock: ThreadLock | ProcessLock | nullcontext | None = None,
        continue_on_failed_credentials: bool = True,
    ):
        super().__init__(
            igvf_mode=igvf_mode,
            submission=submission,
            dry_run=dry_run,
            no_log_file=True,
        )

        self.session = requests.Session()
        self.session.auth = self.auth
        self.logger = ParallelLogger.new(logger=self.debug_logger, lock=lock)
        self.continue_on_failed_credentials = continue_on_failed_credentials

    @classmethod
    def new(
        cls,
        igvf_mode: IgvfMode,
        submission=False,
        dry_run=False,
        lock: ThreadLock | ProcessLock | nullcontext | None = None,
        continue_on_failed_credentials: bool = True,
    ) -> PConnection:
        connection = cls(
            igvf_mode=igvf_mode,
            submission=submission,
            dry_run=dry_run,
            continue_on_failed_credentials=continue_on_failed_credentials,
        )
        if isinstance(lock, ProcessLock):
            # in a process-parallel environment, may need to clean up the
            # root logger again
            utils.fix_igvf_logging()

        return connection

    def set_submission(self, status: bool):
        self.submission = status

    def check_dry_run(self) -> bool:
        return self.dry_run

    @classmethod
    def _props_equal[T](cls, prop1: T, prop2: T) -> bool:
        if prop1 == prop2:
            return True
        else:
            match prop1, prop2:
                case str(p1), str(p2):
                    return cls._props_equal(p1.split(","), p2.split(","))
                case str(p1), list(p2):
                    return cls._props_equal(p1.split(","), p2)
                case list(p1), str(p2):
                    return cls._props_equal(p1, p2.split(","))
                case list(p1), list(p2):
                    return sorted(p1) == sorted(p2)
                case _:
                    return False

    def _patch_in_post(
        self,
        record_id: str,
        profile: IgvfSchema,
        payload: dict[str, object],
        existing_record: dict[str, object],
    ) -> None:
        changed = False
        changed_props = set()
        added_props = set()
        removed_props = set()
        for key, value in payload.items():
            if key in {self.IGVFID_KEY, self.PROFILE_KEY}:
                continue
            if existing_record.get(key, None) is None:
                changed = True
                added_props.add(key)
                existing_record[key] = value
            elif not self._props_equal(existing_record[key], value):
                changed = True
                existing_record[key] = value
                changed_props.add(key)
        to_remove_props = set(existing_record.keys()).difference(payload.keys())
        for key in to_remove_props:
            prop = profile.get_property_from_name(key)
            if prop.is_required or prop.is_not_submittable or prop.is_read_only:
                continue  # we don't remove this property
            # Then it is safe to remove this property.
            changed = True
            existing_record.pop(key)
            removed_props.add(key)

        if not changed:
            self.logger.info(f"No changes to record for '{record_id}', leaving as-is.")
            return

        url = iuu.url_join([self.igvf_mode.url, record_id.lstrip("/")])
        response = self.session.put(
            url,
            timeout=iu.TIMEOUT,
            headers=iuu.REQUEST_HEADERS_JSON,
            json=existing_record,
            verify=False,
        )
        response_json = response.json()

        if response.ok:
            self.logger.info(
                f"Successfully PUT {record_id}.\n"
                f"\tadded {','.join(added_props)}; changed {','.join(changed_props)}; removed {','.join(removed_props)}\n"
            )
        else:
            self.logger.error(
                f"Failed to PUT {record_id}\n"
                f"\tadded {','.join(added_props)}; changed {','.join(changed_props)}; removed {','.join(removed_props)}\n"
                f"{iuu.print_format_dict(response_json)}"
            )
            response.raise_for_status()

    def post(
        self,
        payload: dict[str, object],
        require_aliases: bool = True,
        upload_file: bool = True,
        return_original_status_code: bool = False,
        truncate_long_strings_in_payload_log: bool = False,
        upload_duplicate: bool = True,
    ):
        """POST a record to the Portal.

        Requires that you include in the payload the non-schematic key ``self.PROFILE_KEY`` to
        designate the name of the IGVF object profile that you are submitting to, or the
        actual `@id` property itself.

        If the `lab` property isn't present in the payload, then the default will be set to the value
        of the `IGVF_LAB` environment variable. Similarly, if the `award` property isn't present, then the
        default will be set to the value of the `IGVF_AWARD` environment variable.

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

        Returns:
            `dict`: The JSON response from the POST operation, or the existing record if it already
            exists on the Portal (where a GET on any of it's aliases, when provided in the payload,
            finds the existing record). If `return_original_status_code=True`, then will
            return a `tuple` of the above `dict` and an `int` corresponding to the
            status code on POST of the initial payload.

        Raises:
            igvf_utils.exceptions.AwardPropertyMissing: The `award` property isn't present in the payload and there isn't a
                default set by the environment variable `IGVF_AWARD`.
            igvf_utils.exceptions.LabPropertyMissing: The `lab` property isn't present in the payload and there isn't a
                default set by the environment variable `IGVF_LAB`.
            igvf_utils.exceptions.MissingAlias: The argument 'require_aliases' is set to True and
                the 'aliases' property is missing in the payload or is empty.
            requests.exceptions.HTTPError: The return status is not ok.

        Side effects:
            self.PROFILE_KEY will be popped out of the payload if present, otherwise, the key "@id"
            will be popped out. Furthermore, self.IGVFID_KEY will be popped out if present in the payload.
        """
        self.logger.debug("\nIN post().")
        # Make sure we have a payload that can be converted to valid JSON, and
        # tuples become arrays, ...
        payload = json.loads(json.dumps(payload))
        profile = self.get_profile_from_payload(payload)
        payload[self.PROFILE_KEY] = profile.name
        url = iuu.url_join([self.igvf_mode.url, profile.name])
        if self.IGVFID_KEY in payload:
            # Shouldn't be here, unless maybe a PATCH was attempted and the record didn't exist, so
            # a POST was then attempted.
            payload.pop(self.IGVFID_KEY)
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
        try:
            payload.pop(self.PROFILE_KEY)
        except KeyError:
            pass
        try:
            payload.pop("@id")
        except KeyError:
            pass

        no_alias = False  # Use this to check later if doing a GET
        aliases = payload.get(iu.ALIAS_PROP_NAME)
        if not aliases:
            if not profile.has_alias or not require_aliases:
                aliases = ["N/A"]
                no_alias = True
            else:
                raise MissingAlias(
                    (
                        "Missing property '{}' in payload {}. This is required by default for the profiles"
                        " that include this property, and can be disabled by setting the `require_aliases`"
                        " argument to False in the call to this method, being `igvf_utils.connection.Connection.post()`."
                        " If using the iu_register.py script, this can also be disabled by passing the --no-aliases option."
                    ).format(iu.ALIAS_PROP_NAME, payload)
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
            self.logger.error("Payload is: {}".format(iuu.print_format_dict(payload)))
            self.logger.error(validation_error[0])  # The top-level validation message
            if validation_error[1]:  # The validation context can be empty
                self.logger.error(iuu.print_format_dict(validation_error[1]))
            raise Exception(iuu.print_format_dict(validation_error[0]))

        self.logger.debug(
            (
                f"<<<<<< POST {profile.name} record {aliases[0]} To IGVF database with URL {url} and this payload:\n"
                f"{iuu.print_format_dict(payload, truncate_long_strings=truncate_long_strings_in_payload_log)}"
            )
        )

        if self.check_dry_run():
            return {}
        response = self.session.post(
            url,
            timeout=iu.TIMEOUT,
            headers=iuu.REQUEST_HEADERS_JSON,
            json=payload,
            verify=False,
        )
        response_json = response.json()
        original_status_code = response.status_code

        if response.ok:
            self.logger.debug("Success.")
            response_json = response_json["@graph"][0]
            # Some objects don't have an accession, i.e. replicates.
            # in original code "record_id" was frequently called "encid" (presumably encode id?)
            record_id: str = response_json.get("accession", response_json["uuid"])
            self.logger.debug(f"Object posted with identifier: {record_id}")
            self._log_post(aliases=aliases, dacc_id=record_id)
            # Run 'after' hooks:
            self.after_submit_hooks(
                record_id, profile.name, method=self.POST, upload_file=upload_file
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
            #   'detail': "Keys conflict: [('file:paired_with', 'f39320d9-0970-4369-b680-5965a5e85b6f')]",
            #   'description': 'There was a conflict when trying to complete your request.',
            #   'code': 409,
            #   '@type': ['HTTPConflict', 'Error'],
            #   'title': 'Conflict',
            #   'status': 'error'}
            # }
            #
            if no_alias:
                self.logger.warning(response_json)
                response.raise_for_status()
            else:
                existing_record = self.get(
                    rec_ids=aliases, ignore404=False, frame="edit"
                )
                if not existing_record:
                    self.logger.warning(response_json)
                    response.raise_for_status()
                else:
                    if upload_duplicate:
                        record_id: str = existing_record.get(
                            "@id", existing_record.get("accession", aliases[0])
                        )
                        self.logger.warning(
                            f"Conflict when POSTing {record_id}, will patch any differences."
                        )
                        self._patch_in_post(
                            record_id=record_id,
                            profile=profile,
                            payload=payload,
                            existing_record=existing_record,
                        )
                        self.after_submit_hooks(
                            record_id,
                            profile.name,
                            method=self.POST,
                            upload_file=upload_file,
                        )
                    else:
                        self.logger.error(
                            f"Will not POST '{aliases[0]}' since it already exists with aliases '{existing_record['aliases']}'."
                        )
                    if return_original_status_code is True:
                        return (existing_record, original_status_code)
                    return existing_record

        else:
            self.logger.error(
                f"Failed to POST {aliases[0]}\n{iuu.print_format_dict(response_json)}"
            )
            response.raise_for_status()

    def _already_uploaded(
        self,
        upload_url: str,
        file_path: Path,
        file_id: str,
        aws_creds: dict | None = None,
    ) -> bool:
        if aws_creds is None:
            upload_credentials = self.get_upload_credentials(file_id)
            aws_creds = self.extract_aws_upload_credentials(upload_credentials)
        bucket, key = utils.parse_s3_uri(upload_url)
        # TODO: replace with boto3 call
        result = subprocess.run(
            f"aws s3api head-object --bucket '{bucket}' --key '{key}'",
            shell=True,
            capture_output=True,
            env={**os.environ, **aws_creds},
        )

        if result.returncode != 0:
            self.logger.info(f"Upload url '{upload_url}' doesn't already exist.")
            return False
        else:
            result = json.loads(result.stdout)
            remote_md5_sum = result.get("Metadata", {}).get("md5sum", "")
            local_md5_sum = utils.md5sum(file_path)
            # this is unneccessary: interrupted/failed uploads do not result in remote obects,
            # so the only way the metadata could be present is if the object uploads successfully
            # remote_crc64_nvme_checksum = result.get("ChecksumCRC64NVME", "")
            # local_crc64_nvme_checksum = utils.crc64_nvme_checksum(file_path)
            alread_exists = remote_md5_sum == local_md5_sum
            if alread_exists:
                self.logger.info(
                    f"Upload url '{upload_url}' md5_sum matches local file, will not re-upload."
                )
            else:
                self.logger.info(
                    f"Upload url '{upload_url}' md5_sum does NOT match local file, will re-upload."
                )
            return alread_exists

    def upload_file(
        self,
        file_id: str,
        file_path: str | Path | None = None,
        set_md5sum: str | bool | None = None,
    ):
        """
        Uploads a file to the Portal for the indicated file record. The file to upload can be
        specified by setting the `file_path` parameter, or by using the value of the IGVF file
        profile's `submitted_file_name` property of the given file object represented by the
        `file_id` parameter. The file to upload can be from any of the following sources:

          1. Path to a local file,
          2. S3 object, or
          3. Google Storage object

        For the AWS option above, the user must set the proper AWS keys, see the `wiki documentation`_.

        If the dry-run feature is enabled, then this method will return prior to launching the
        upload command.

        Args:
            file_id: `str`. An identifier of a `file` record on the IGVF Portal.
            file_path: `str`. The local path to the file to upload, or an S3 object (i.e s3://mybucket/test.txt),
              or a Google Storage object (i.e. gs://mybucket/test.txt).
              If not set, defaults to `None` in which case the local file path will be extracted from the
              record's `submitted_file_name` property.
            set_md5sum: `bool`. True means to also calculate the md5sum and set the file record's `md5sum`
              property on the Portal (this currently is only implemented for local files and S3; not yet GCP).
              This will always take place whenever the property isn't yet set.
              Furthermore, setting to True will also cause the `file_size` property to be set.
              Normally these two properties would already be set as they are required in the *file* profile,
              however, if the wrong file was originally uploaded, then they must be reset when
              uploading a new file.

        Raises:
            igvf_utils.exceptions.FileUploadFailed: The return code of the AWS upload command was non-zero.

        .. _`wiki documentation`: https://github.com/IGVF-DACC/igvf_utils/wiki/Configuration#aws-keys
        """
        self.logger.debug("\nIN upload_file()\n")
        # upload_credentials = self.get_upload_credentials(file_id) # Don't use this - they may have expired.
        try:
            upload_credentials = self.regenerate_aws_upload_creds(file_id)
            aws_creds = self.extract_aws_upload_credentials(upload_credentials)
        except requests.exceptions.HTTPError:
            if self.continue_on_failed_credentials:
                self.logger.warning(
                    f"Skipping upload of {file_id} because credentials are not available. It is probably finalized."
                )
                return
            else:
                raise
        file_rec: IgvfRecord = self.get(rec_ids=file_id, ignore404=False)
        match file_path:
            case None:
                try:
                    file_path: Path = Path(
                        file_rec[self.profiles.SUBMITTED_FILE_PROP_NAME]
                    )
                except KeyError:  # submitted_file_name property not set:
                    raise Exception("No file path specified.")
            case str(str_path):
                file_path = Path(str_path)
            case Path():
                pass
            case _:
                raise ValueError(f"Invalid type {type(file_path)} for file_path.")

        upload_url = aws_creds["UPLOAD_URL"]

        if self._already_uploaded(
            upload_url=upload_url,
            file_path=file_path,
            file_id=file_id,
            aws_creds=aws_creds,
        ):
            return

        cmd = f"aws s3 cp {file_path} {upload_url} --metadata 'md5sum={utils.md5sum(file_path)}'"
        self.logger.info(f"Uploading {file_path} to {upload_url}")
        self.logger.debug(f"Running command '{cmd}'.")
        if self.check_dry_run():
            return
        # TODO: replace with boto3 call
        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            env=os.environ.update(aws_creds),
        )
        if result.returncode != 0:
            self.logger.error(
                f"Subprocess command '{cmd}' failed with return code '{result.returncode}'.\n"
                f"Stdout is '{result.stdout.decode('utf-8')}'\n."
                f"Stderr is '{result.stderr.decode('utf-8')}'."
            )
            raise RuntimeError(f"Failed to upload file '{file_path}' for {file_id}.")
        self.logger.info(f"Successfully uploaded {file_path} to {upload_url}.")

    def get(self, rec_ids, database=False, ignore404=True, frame=None):
        """GET a record from the Portal.

        Looks up a record in the Portal and performs a GET request, returning the JSON serialization of
        the object. You supply a list of identifiers for a specific record, and the Portal will be
        searched for each identifier in turn until one is either found or the list is exhausted.

        Args:
            rec_ids: `str` or `list`. Must be a `list` if you want to supply more than one identifier.
                For a few example identifiers, you can use a uuid, accession, ..., or even the value of
                a record's `@id` property.
            database: `bool`. If True, then search the database directly instead of the Elasticsearch.
                 indices. Always True when in submission mode (`self.submission` is True).
            frame: `str`. A value for the frame query parameter, i.e. 'object', 'edit'.
            ignore404: `bool`. Only matters when none of the passed in record IDs were found on the
                Portal.  In this case, If set to `True`, then None will be returned.
                If set to `False`, then an Exception will be raised.


        Returns:
            `dict`: The JSON response. Will be empty if no record was found AND ``ignore404=True``.

        Raises:
            `Exception`: If the server responds with a FORBIDDEN status.
            `requests.exceptions.HTTPError`: The status code is not ok, and the
                cause isn't due to a 404 (not found) status code when ``ignore404=True``.
        """
        if self.submission:
            database = True
        if isinstance(rec_ids, str):
            rec_ids = [rec_ids]
        status_codes = {}  # key is return code, value is the record ID
        for r in rec_ids:
            r = r.strip("/")
            url = iuu.url_join([self.igvf_mode.url, r, "?format=json"])
            if database:
                url += "&datastore=database"
            if frame:
                url += "&frame={frame}".format(frame=frame)
            self.logger.debug(f">>>>>>GET {r} From DACC with URL {url}")
            response = self.session.get(
                url,
                timeout=iu.TIMEOUT,
                headers=iuu.REQUEST_HEADERS_JSON,
                verify=False,
            )
            if response.ok:
                return response.json()
            status_codes[response.status_code] = r

        if requests.codes.FORBIDDEN in status_codes:
            raise Exception(
                "Access to IGVF record {} is forbidden".format(
                    status_codes[requests.codes.FORBIDDEN]
                )
            )
        elif requests.codes.NOT_FOUND in status_codes:
            self.logger.debug("NOT FOUND")
            if ignore404:
                return None
        # At this point in the code, the response is not okay.
        # Raise the error for last response we got:
        response.raise_for_status()
