import json

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

from igvf_portal import utils
from igvf_portal.enums import IgvfMode


class PConnection(Connection):
    def __init__(
        self,
        igvf_mode: IgvfMode,
        submission: bool = False,
        dry_run: bool = False,
        log_file: bool = False,
    ):
        super().__init__(
            igvf_mode=igvf_mode,
            submission=submission,
            dry_run=dry_run,
            no_log_file=not log_file,
        )

    @classmethod
    def new(cls, igvf_mode: IgvfMode, submission=False, dry_run=False) -> PConnection:
        connection = cls(igvf_mode=igvf_mode, submission=submission, dry_run=dry_run)
        utils.fix_igvf_logging(debug_logger=connection.debug_logger)
        return connection

    def set_submission(self, status: bool):
        self.submission = status

    def check_dry_run(self) -> bool:
        return self.dry_run

    def post(
        self,
        payload: object,
        require_aliases: bool = True,
        upload_file: bool = True,
        return_original_status_code: bool = False,
        truncate_long_strings_in_payload_log: bool = False,
        upload_duplicate: bool = False,
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
        self.debug_logger.debug("\nIN post().")
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
        self.debug_logger.debug("Validating the payload against the schema")
        validation_error = iuu.err_context(
            payload=payload,
            schema=self.profiles.get_profile_from_id(profile.name).schema,
        )
        if validation_error:
            self.debug_logger.error(
                f"Invalid schema instance of the {profile.name} profile."
            )
            self.debug_logger.error(
                "Payload is: {}".format(iuu.print_format_dict(payload))
            )
            self.debug_logger.error(
                validation_error[0]
            )  # The top-level validation message
            if validation_error[1]:  # The validation context can be empty
                self.debug_logger.error(iuu.print_format_dict(validation_error[1]))
            raise Exception(iuu.print_format_dict(validation_error[0]))

        self.debug_logger.debug(
            (
                f"<<<<<< POST {profile.name} record {aliases[0]} To IGVF database with URL {url} and this payload:\n"
                f"{iuu.print_format_dict(payload, truncate_long_strings=truncate_long_strings_in_payload_log)}"
            )
        )

        if self.check_dry_run():
            return {}
        response = requests.post(
            url,
            auth=self.auth,
            timeout=iu.TIMEOUT,
            headers=iuu.REQUEST_HEADERS_JSON,
            json=payload,
            verify=False,
        )
        # response_json = response.json()["@graph"][0]
        response_json = response.json()
        original_status_code = response.status_code

        if response.ok:
            self.debug_logger.debug("Success.")
            response_json = response_json["@graph"][0]
            # Some objects don't have an accession, i.e. replicates.
            encid: str = response_json.get("accession", response_json["uuid"])
            self.debug_logger.debug(f"Object posted with identifier: {encid}")
            self._log_post(aliases=aliases, dacc_id=encid)
            # Run 'after' hooks:
            self.after_submit_hooks(
                encid, profile.name, method=self.POST, upload_file=upload_file
            )
            if return_original_status_code is True:
                return (response_json, original_status_code)
            return response_json
        elif response.status_code == requests.codes.CONFLICT:
            self.debug_logger.error(response_json)
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
                response.raise_for_status()
            else:
                existing_record = self.get(rec_ids=aliases, ignore404=True)
                if not existing_record:
                    response.raise_for_status()
                else:
                    self.debug_logger.error(
                        f"Will not POST '{aliases[0]}' since it already exists with aliases '{existing_record['aliases']}'."
                    )
                    if upload_duplicate:
                        encid: str = existing_record.get(
                            "accession", existing_record["uuid"]
                        )
                        self.after_submit_hooks(
                            encid,
                            profile.name,
                            method=self.POST,
                            upload_file=upload_file,
                        )
                    if return_original_status_code is True:
                        return (existing_record, original_status_code)
                    return existing_record

        else:
            self.debug_logger.error(f"Failed to POST {aliases[0]}")
            self.debug_logger.error("<<<<<< DACC POST RESPONSE: ")
            self.debug_logger.error(iuu.print_format_dict(response_json))
            response.raise_for_status()
