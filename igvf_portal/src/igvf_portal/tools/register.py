import csv
import json
import multiprocessing
import re
import time
from collections.abc import Iterable, Iterator
from concurrent.futures import ProcessPoolExecutor
from functools import wraps  # for retry function
from pathlib import Path
from typing import Final, cast

import igvf_utils.utils as iuu
from igvf_utils.profiles import IgvfSchema

from igvf_portal import utils
from igvf_portal.connection import PConnection
from igvf_portal.constants import VERSION
from igvf_portal.enums import IgvfMode, LogLevel
from igvf_portal.register_config import RegisterConfig

_STR_REGEX: Final[re.Pattern] = re.compile(r'\'|"')

_CONNECTION: PConnection
_REGISTER_CONFIG: RegisterConfig

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###
# © 2018 The Board of Trustees of the Leland Stanford Junior University
# Nathaniel Watson
# nathankw@stanford.edu
###

"""
Given a tab-delimited, JSON, or JSONL input file containing one or more records belonging to one
of the profiles listed on the IGVF Portal (such as https://sandbox.igvf.org/profiles/document.json),
either POSTS or PATCHES the records. The default is to POST each record; to PATCH instead, see
the ``--patch`` option.

When POSTING file records, the md5sum of each file will be calculated for you if you haven't
already provided the `md5sum` property. Then, after the POST operation completes, the actual file
will be uploaded to AWS S3. In order for this to work, you must set the `submitted_file_name`
property to the full, local path to your file to upload. Alternatively, you can set
`submitted_file_name` to and existing S3 object, i.e. s3://mybucket/reads.fastq.

Note that there is a special 'trick' defined in the ``igvf_utils.connection.Connection()``
class that can be taken advantage of to simplify submission under certain profiles.
It concerns the `attachment` property in any profile that employs it, such as the `document`
profile.  The trick works as follows: instead of constructing the `attachment` propery object
value as defined in the schema, simply use a single-key object of the following format::

  {"path": "/path/to/myfile"}

and the `attachment` object will be constructed for you.

|
"""


RECORD_ID_FIELD = "record_id"
"""RECORD_ID_FIELD is a special field that won't be skipped in the create_payload() function.
It is used when patching objects to indicate the identifier of the record to patch.
"""


# def get_parser():
#     parser = argparse.ArgumentParser(
#         description=__doc__,
#         parents=[igvf_login_parser],
#         formatter_class=argparse.RawTextHelpFormatter,
#     )

#     parser.add_argument(
#         "-d",
#         "--dry-run",
#         action="store_true",
#         help="""
#     Set this option to enable the dry-run feature, such that no modifications are performed on the
#     IGVF Portal.  This is useful if you'd like to inspect the logs or ensure the validity of
#     your input file.""",
#     )

#     parser.add_argument(
#         "--no-aliases",
#         action="store_true",
#         help="""
#     Setting this option is NOT advised. Set this option for doing a POST when your input file
#     doesn't contain an 'aliases' column, even though this property is supported in the corresponding
#     IGVF profile.
#     When POSTING a record to a profile that includes the 'aliases' property, this package requires
#     the 'aliases' property be used for traceability purposes and because without this property,
#     it'll be very easy to create duplicate objects on the Portal.  For example, you can easily
#     create the same biosample as many times as you want on the Portal when not providing an alias.""",
#     )

#     parser.add_argument(
#         "--no-upload-file",
#         action="store_true",
#         help="Don't upload files when POSTing file objects",
#     )

#     parser.add_argument(
#         "-p",
#         "--profile_id",
#         required=True,
#         help="""
#     The ID of the profile to submit to, i.e. use 'document' for
#     https://sandbox.igvf.org/profiles/document.json. The profile will be pulled down for
#     type-checking in order to type-cast any values in the input file to the proper type (i.e. some
#     values need to be submitted as integers, not strings).""",
#     )

#     parser.add_argument(
#         "-i",
#         "--infile",
#         required=True,
#         help="""
#     The JSON, JSONL, or tab-delimited input file.

#     **The tab-delimited file format:**
#     Must have a field-header line as the first line.
#     Any lines after the header line that start with a '#' will be skipped, as well as any empty lines.
#     The field names must be exactly equal to the corresponding property names in the corresponding
#     profile. Non-scematic fields are allowed as long as they begin with a '#'; they will be
#     skipped. If a property has an array data type (as indicated in the profile's documentation
#     on the Portal), the array literals '[' and ']' are optional. Values within the array must
#     be comma-delimited. For example, if a property takes an array of strings, then you can use
#     either of these as the value:

#     1) str1,str2,str3
#     2) [str1,str2,str3]

#     On the other hand, if a property takes a JSON object as a value, then the value you enter must be
#     valid JSON. This is true anytime you have to specify a JSON object. Thus, if you are submitting a
#     genetic_modification and you have two 'introduced_tags' to provide, you can supply them in either
#     of the following two ways:

#     1) {"name": "eGFP", "location": "C-terminal"},{"name": "FLAG","C-terminal"}
#     2) [{"name": "eGFP", "location": "C-terminal"},{"name": "FLAG","C-terminal"}]

#     **The JSON input file**
#     Can be a single JSON object, or an array of JSON objects. Key names must match property names of
#     an IGVF record type (profile).

#     **The JSONL input file**
#     This format is largely similar to the JSON input file, but each row in the file is a single
#     JSON object.

#     **The following applies to all input file formats**
#     When patching objects, you must specify the 'record_id' field to indicate the identifier of the record.
#     Note that this a special field that is not present in the IGVF schema, and doesn't use the '#'
#     prefix to mark it as non-schematic. Here you can specify any valid record identifier
#     (i.e. UUID, accession, alias).

#     Some profiles (most) require specification of the 'award' and 'lab' attributes. These may be set
#     as fields in the input file, or can be left out, in which case the default values for these
#     attributes will be pulled from the environment variables IGVF_AWARD and IGVF_LAB, respectively.
#     """,
#     )

#     parser.add_argument(
#         "-w",
#         "--overwrite-array-values",
#         action="store_true",
#         help="""
#     Only has meaning in combination with the --patch option. When this is specified, it means that
#     any keys with array values will be overwritten on the IGVF Portal with the corresponding value
#     to patch. The default action is to extend the array value with the patch value and then to remove
#     any duplicates.""",
#     )

#     parser.add_argument(
#         "-r",
#         "--remove-property",
#         help="""
#     Only has meaning in combination with the --rm-patch option. Properties specified in this argument
#     will be popped from the record fetched from the IGVF portal. Can specify as comma delimited
#     string.""",
#     )

#     group = parser.add_mutually_exclusive_group()

#     group.add_argument(
#         "--patch",
#         action="store_true",
#         help="""
#     Presence of this option indicates to PATCH an existing DCC record rather than register a new one.""",
#     )

#     group.add_argument(
#         "--rm-patch",
#         action="store_true",
#         help="""
#     Presence of this option indicates to remove a property, as specified by the -r argument,
#     from an existing DACC record, and then PATCH it with the payload specified in -i.""",
#     )

#     parser.add_argument(
#         "--tries",
#         type=int,
#         default=1,
#         help="""
#     Number of times to try before giving up to prevent time out error when doing post or patch on a large set.""",
#     )

#     parser.add_argument(
#         "--delay",
#         type=int,
#         default=5,
#         help="""
#     Initial delay between retries in seconds.""",
#     )

#     parser.add_argument(
#         "--backoff",
#         type=int,
#         default=2,
#         help="""
#     Backoff multiplier, by default will double the delay each retry.""",
#     )
#     return parser


##decorator for preventing time out##
def retry(tries=1, delay=5, backoff=2):
    """Retry calling the decorated function using an exponential backoff.

    http://www.saltycrane.com/blog/2009/11/trying-out-retry-decorator-python/
    original from: http://wiki.python.org/moin/PythonDecoratorLibrary#Retry

    args:
        tries (int): number of times to try (not retry) before giving up
        delay (int): initial delay between retries in seconds
        backoff (int): backoff multiplier e.g. value of 2 will double the delay each retry
    """

    def deco_retry(f):

        @wraps(f)
        def f_retry(*args, **kwargs):
            mtries, mdelay = tries, delay
            while mtries > 1:
                try:
                    return f(*args, **kwargs)
                except Exception as e:
                    print(str(e))
                    msg = "Retrying in %d seconds..." % (mdelay)
                    print(msg)
                    time.sleep(mdelay)
                    mtries -= 1
                    mdelay *= backoff
            return f(*args, **kwargs)

        return f_retry  # true decorator

    return deco_retry


def _init_worker(
    register_config: RegisterConfig,
    schema: IgvfSchema,
) -> None:
    """Initialize _REGISTER_CONFIG and _CONNECTION in worker process.

    Set the profile to the updated schema (so that properties are loaded).
    """
    global _REGISTER_CONFIG
    _REGISTER_CONFIG = register_config
    global _CONNECTION
    _CONNECTION = _REGISTER_CONFIG.connection
    _CONNECTION.profiles.get_profile_from_id(register_config.profile_id)
    _CONNECTION.profiles._profiles[register_config.cleaned_profile_id] = schema


def _check_valid_json(prop: str, val: str, line_number: int) -> object:
    """
    Runs json.loads(val) to ensure valid JSON.

    Args:
        val: str. A string load as JSON.
        prop: str. Name of the schema property/field that stores the passed in val.
        row_count: int. The line number from the input file that is currently being processed.

    Raises:
        ValueError: The input is malformed JSON.
    """

    # Don't try to break down the individual pieces of a nested object. That will be too complex for this script, and will also
    # be too complex for the end user to try and represent in some flattened way. Thus, require the end user to supply proper JSON
    # for a nested object.
    json_val = json.loads(val)
    if isinstance(json_val, list):
        for item in json_val:
            if not isinstance(item, dict):
                raise ValueError(
                    f"Error: Invalid JSON in field '{prop}', row '{line_number}'"
                )
    return json_val


def _typecast(field_name, value, data_type, line_num):
    """
    Converts the value to the specified data type. Used to convert string representations of integers
    in the input file to integers, and string representations of booleans to booleans.

    Args:
        field_name: The name of the field in the input file whose value is being potentially typecast.
            Used only in error messages.
        value: The value to potentially typecast.
        data_type: Specifies the data type of field_name as indicated in the IGVF profile.
        line_num: The current line number in the input file. Used only in error messages.
    """
    match data_type:
        case "integer":
            return int(value)
        case "number":
            # JSON Schema says that a number can by any numeric type.
            # First check if integer, if not, treat as float.
            try:
                return int(value)
            except ValueError:
                # This will be raised if trying to convert a string representation of a float to an int.
                return float(value)
        case "boolean":
            value = value.lower()
            if value not in ["true", "false"]:
                raise ValueError(
                    f"Can't convert value '{value}' in field '{field_name}' on line {line_num} to data type '{data_type}'."
                )
            return value == "true"
        case _:
            return value


def _get_field_val(
    field: str, val: str, schema: IgvfSchema, line_number: int
) -> object | None:
    """Get value for payload field from input str representation."""
    if len(val) == 0:
        return None
    if field == RECORD_ID_FIELD:
        return val
    field_schema = schema.get_property_from_name(field).schema
    match field_schema["type"]:
        case "object":
            # Must be proper JSON
            return _check_valid_json(field, val, line_number)
        case "array":
            item_val_type = field_schema["items"]["type"]
            if item_val_type == "object":
                # Must be valid JSON
                # Check if user supplied optional JSON array literal. If not, I'll add it.
                if not val.startswith("["):
                    val = "[" + val
                if not val.endswith("]"):
                    val += "]"
                return _check_valid_json(field, val, line_number)
            else:
                # User is allowed to enter values in string literals. I'll remove them if I find them,
                # since I'm splitting on the ',' to create a list of strings anyway:
                val = _STR_REGEX.sub("", val)
                # Remove optional JSON array literal since I'm tokenizing and then converting
                # to an array regardless.
                if val.startswith("["):
                    val = val[1:]
                if val.endswith("]"):
                    val = val[:-1]
                val_list = [x.strip() for x in val.split(",")]
                # Type cast tokens if need be, i.e. to integers:
                return [
                    _typecast(
                        field_name=field,
                        value=x,
                        data_type=item_val_type,
                        line_num=line_number,
                    )
                    for x in val_list
                    if x
                ]
        case _ as schema_val_type:  # schema type not object or array
            return _typecast(
                field_name=field,
                value=val,
                data_type=schema_val_type,
                line_num=line_number,
            )


def _iter_payloads_from_tsv(
    schema: IgvfSchema, infile: Path
) -> Iterator[dict[str, object]]:
    """
    Generates the payload for each row in 'infile'.

    Args:
        schema: IgvfSchema. The schema of the objects to be submitted.
        infile - str. Path to input file.

    Yields  : dict. The payload that can be used to either register or patch the metadata for each row.
    """
    # Fetch the schema from the IGVF Portal so we can set attr values to the
    # right type when generating the payload (dict).
    schema_props: list[str] = [prop.name for prop in schema.properties]
    with infile.open("rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Unable to parse header of {infile}")
        field_names: set[str] = set()
        for field in reader.fieldnames:
            if field.startswith("#"):  # non-schema field, don't use it
                continue
            if field not in schema_props:
                if field != RECORD_ID_FIELD:
                    raise ValueError(
                        f"Unknown field name '{field}', which is not registered as a property in the specified schema at {schema.name}."
                    )
            field_names.add(field)

        for line_number, line in enumerate(reader, start=2):
            payload = {PConnection.PROFILE_KEY: schema.name}
            for field_name in field_names:
                val = _get_field_val(
                    field=field_name,
                    val=line[field_name].strip(),
                    line_number=line_number,
                    schema=schema,
                )
                if val is not None:
                    payload[field_name] = val

            yield payload


def _remove_and_patch_payload(payload: dict[str, object]) -> str:
    """Execute the remove_and_patch function for the payload, using options from _REGISTER_CONFIG."""
    record_id = payload.get(RECORD_ID_FIELD, False)
    if not record_id:
        raise ValueError(
            "Can't patch payload {} since there isn't a '{}' field indicating an identifier for the record to be PATCHED.".format(
                iuu.print_format_dict(payload), RECORD_ID_FIELD
            )
        )
    payload.pop(RECORD_ID_FIELD)
    payload.update({_CONNECTION.IGVFID_KEY: record_id})

    _REGISTER_CONFIG.retry()(_CONNECTION.remove_and_patch)(
        props=_REGISTER_CONFIG.remove_properties,
        patch=payload,
        extend_array_values=not _REGISTER_CONFIG.overwrite_array_values,
    )

    return cast(list[str], payload["aliases"])[0]


def _patch_payload(payload: dict[str, object]) -> str:
    """Execute the patch function for the payload, using options from _REGISTER_CONFIG."""
    record_id = payload.get(RECORD_ID_FIELD, False)
    if not record_id:
        raise ValueError(
            "Can't patch payload {} since there isn't a '{}' field indicating an identifier for the record to be PATCHED.".format(
                iuu.print_format_dict(payload), RECORD_ID_FIELD
            )
        )
    payload.pop(RECORD_ID_FIELD)
    payload.update({_CONNECTION.IGVFID_KEY: record_id})

    _REGISTER_CONFIG.retry()(_CONNECTION.patch)(
        payload=payload, extend_array_values=not _REGISTER_CONFIG.overwrite_array_values
    )

    return cast(list[str], payload["aliases"])[0]


def _post_payload(payload: dict[str, object]) -> str:
    """Execute the post function for the payload, using options from _REGISTER_CONFIG."""

    _REGISTER_CONFIG.retry()(_CONNECTION.post)(
        payload,
        require_aliases=True,
        upload_file=_REGISTER_CONFIG.upload_file,
        upload_duplicate=_REGISTER_CONFIG.upload_duplicate,
    )

    return cast(list[str], payload["aliases"])[0]


def register(
    *,
    infile: Path,
    profile_id: str,
    dry_run: bool = False,
    igvf_mode: IgvfMode = IgvfMode.prod,
    patch: bool = False,
    overwrite_array_values: bool = False,
    remove_property: Iterable[str] = (),
    tries: int = 2,
    delay: float = 5.0,
    backoff: float = 2.0,
    upload_file: bool = True,
    upload_duplicate: bool = False,
    num_workers: int = 12,
    update_secs: float = 30.0,
    log_level: LogLevel = LogLevel.info,
):
    """Register data with the IGVF Portal.

    Args:
        infile: The JSON, JSONL, or tab-delimited input file.
        profile_id: @id of data to register
        dry_run: if True, don't actually modify the portal or upload data
        igvf_mode: Which IGVF server to use: "prod", "staging", or "sandbox"
        patch: if True, patch data instead of POSTing new data
        overwrite_array_values: If True, when patching data with array values, overwrite old values
          with new valus. If False, extend old data with new values.
        remove_property: Any properties specified here will be removed from the record before patching.
        tries: Number of times to try (not retry) in case of failure.
        delay: Wait time in seconds before retrying after failure.
        backoff: Multipliciative scale to delay for each successive failure.
        upload_file: If True, actually upload files when posting.
        upload_duplicate: If True, upload file even if a duplicate record exists. Used for uploading files
            when a record POSTed but the upload failed.
        num_workers: Number of connections to the IGVF Portal. If <= 0, use one per CPU.
        update_secs: Time in seconds between progress logs
        log_level: Log level for output
    """
    utils.check_access_keys()
    utils.fix_igvf_logging()
    logger = utils.get_logger_from_file(__file__, level=log_level.value)
    logger.info(f"Version: {VERSION}")

    register_config = RegisterConfig(
        igvf_mode=igvf_mode,
        dry_run=dry_run,
        profile_id=profile_id,
        num_tries=tries,
        delay=delay,
        backoff=backoff,
        overwrite_array_values=overwrite_array_values,
        remove_properties=[
            p for property in remove_property for p in property.split(",")
        ],
        upload_file=upload_file,
        upload_duplicate=upload_duplicate,
    )

    # Get connection into submit mode:
    connection = register_config.connection
    schema: IgvfSchema = connection.profiles.get_profile_from_id(profile_id)

    register_func = (
        _remove_and_patch_payload
        if register_config.rm_patch
        else _patch_payload
        if patch
        else _post_payload
    )

    payloads = list(_iter_payloads_from_tsv(schema=schema, infile=infile))

    if num_workers <= 0:
        num_workers = multiprocessing.cpu_count()
    with ProcessPoolExecutor(
        max_workers=num_workers,
        initializer=_init_worker,
        initargs=(register_config, schema),
    ) as executor:
        mapped_futures = executor.map(
            register_func, payloads, buffersize=2 * num_workers
        )
        num_futures = len(payloads)
        t_start = time.time()
        t_last = t_start
        for idx, record_id in enumerate(mapped_futures, start=1):
            logger.debug(f"Processed {record_id}")
            t_now = time.time()
            if t_now - t_last > update_secs:
                t_last = t_now
                logger.info(f"Processed {idx} / {num_futures} records.")

    duration = time.time() - t_start
    logger.info(f"Processed {num_futures} records in {duration:.1f} sec.")
