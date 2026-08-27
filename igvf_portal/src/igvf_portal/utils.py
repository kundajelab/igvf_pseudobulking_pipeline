import csv
import gzip
import hashlib
import logging
import os
import re
import sys
import time
import unicodedata
from collections.abc import (
    Collection,
    Generator,
    Iterable,
    Iterator,
)
from contextlib import contextmanager
from io import TextIOWrapper
from pathlib import Path
from typing import Callable, Literal, TextIO

import igvf_utils
import requests
from igvf_client import ApiClient, Configuration, IgvfApi

from igvf_portal.enums import IgvfMode
from igvf_portal.types import FromTypedDict, IgvfRecord


def setup_logger(logger: logging.Logger, level: int = logging.INFO) -> None:
    logger.handlers.clear()
    logger.setLevel(level)


def get_logger(name: str, level: int = logging.INFO) -> logging.Logger:
    logger = logging.getLogger(name)
    setup_logger(logger, level)
    return logger


def get_logger_from_file(file_name: str, level: int = logging.INFO) -> logging.Logger:
    _, package_name, _, tool_name = file_name.split(".", 1)[0].rsplit("/", 3)
    logger = logging.getLogger(
        f"{package_name.replace('_', '-')} {tool_name.replace('_', '-')}"
    )
    setup_logger(logger, level)
    return logger


def fix_igvf_logging(
    level: int = logging.INFO, debug_logger: logging.Logger | None = None
):
    root = logging.getLogger()
    setup_logger(root, level=level)
    ch = logging.StreamHandler(stream=sys.stderr)
    f_formatter = logging.Formatter(
        "%(asctime)s %(name)s %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    ch.setLevel(level)
    ch.setFormatter(f_formatter)
    root.addHandler(ch)
    if debug_logger is None:
        debug_logger = logging.getLogger("iu_debug")
    setup_logger(debug_logger, level=level)
    error_logger = logging.getLogger("iu_error")
    error_logger.propagate = False
    error_logger.handlers.clear()
    error_logger.addHandler(logging.NullHandler())


def check_access_keys():
    if "IGVF_API_KEY" not in os.environ or "IGVF_SECRET_KEY" not in os.environ:
        raise ValueError("IGVF_API_KEY and IGVF_SECRET_KEY must be set in environment")


def open_igvf_api(igvf_mode: IgvfMode | str) -> IgvfApi:
    """Return an IgvfApi object"""
    check_access_keys()
    _igvf_mode = igvf_mode if isinstance(igvf_mode, IgvfMode) else IgvfMode[igvf_mode]
    config = Configuration(
        access_key=os.environ["IGVF_API_KEY"],
        secret_access_key=os.environ["IGVF_SECRET_KEY"],
        host=_igvf_mode.portal_api_url,
    )
    client = ApiClient(config)
    return IgvfApi(client)


def lookup_ontology_by_cl_id(cl_id: str) -> list[dict[str, object]] | None:
    """Retrieve the full OLS record for a term given its CL_id or short_form ID. Example: CL:0011026"""
    response = requests.get(
        "https://www.ebi.ac.uk/ols4/api/terms",
        params={"short_form": cl_id.replace(":", "_")},
    )
    response.raise_for_status()
    data = response.json()

    # The response is a paginated object; terms are in the embedded "terms" list
    return data.get("_embedded", {}).get("terms", None)


@contextmanager
def maybe_gzipped(maybe_gzipped: Path, mode: Literal["w", "r"]) -> Generator[TextIO]:
    """Get a text reader from a path that may or may not be gzipped."""
    if maybe_gzipped.name.endswith(".gz"):
        if mode == "r":
            with gzip.open(f"{maybe_gzipped}", "rb") as f_gzip:
                yield TextIOWrapper(f_gzip)
        else:
            with gzip.open(f"{maybe_gzipped}", "wb") as f_gzip:
                yield TextIOWrapper(f_gzip)
    else:
        yield maybe_gzipped.open("rt" if mode == "r" else "wt")


def get_sep_from_path(csv_or_tsv: Path) -> Literal["\t", ","]:
    """Guess the file separator from the path."""
    match csv_or_tsv.suffixes:
        case [*_parts, ".tsv"] | [*_parts, ".tsv", ".gz"]:
            return "\t"
        case [*_parts, ".csv"] | [*_parts, ".csv", ".gz"]:
            return ","
        case _:
            raise ValueError(
                f"Could not infer separator for file {csv_or_tsv} with suffix {csv_or_tsv.suffix}. Please "
                "specify separator with 'sep' argument."
            )


def iter_csv_rows(
    csv_path, sep: str | None = None, required_columns: Collection[str] | None = None
) -> Iterator[dict[str, str]]:
    if sep is None:
        sep = get_sep_from_path(csv_path)
    with maybe_gzipped(csv_path, mode="r") as f:
        reader = csv.DictReader(f, delimiter=sep)
        if reader.fieldnames is None:
            raise ValueError(f"No fieldnames in {csv_path}")
        if required_columns is not None:
            reader.fieldnames = [name.strip() for name in reader.fieldnames]
            missing = set(required_columns).difference(reader.fieldnames)
            if len(missing) > 0:
                raise ValueError(
                    f"{csv_path} is missing required columns: {','.join(missing)}"
                )
        yield from reader


def iter_pseudobulk_dirs(pseudobulk_dir: Path) -> Iterator[Path]:
    """List all the sub-folders in the pseudobulk_dir. Each should correspond to a unique pseudobulk."""
    for folder in sorted(pseudobulk_dir.iterdir()):
        if folder.is_dir():
            yield folder


def parse_s3_uri(s3_uri: str) -> tuple[str, str]:
    """Get bucket and object name from AWS S3 URI."""
    if not s3_uri.startswith("s3://"):
        raise ValueError(f"'{s3_uri}' is not a valid S3 path.")
    _, _, bucket_name, object_key = s3_uri.split("/", 3)
    return bucket_name, object_key


def get_record_http_url(
    file_record: IgvfRecord,
    igvf_mode: IgvfMode,
    igvf_portal_region: str = "us-west-2",
) -> str:
    """Get a public HTTPS URL for downloading the specified file record."""
    s3_uri = file_record["s3_uri"]
    if s3_uri.startswith("s3://igvf-public/"):
        # it's a public S3 URL, construct equivalent https URL and return it
        _, _, bucket_name, object_key = s3_uri.split("/", 3)
        return (
            f"https://{bucket_name}.s3.{igvf_portal_region}.amazonaws.com/{object_key}"
        )
    else:
        # it's a private S3 URL, but we can request a temporary download via the "href" field
        response = requests.head(
            url=f"{igvf_mode.portal_api_url}{file_record['href']}",
            auth=(os.environ["IGVF_API_KEY"], os.environ["IGVF_SECRET_KEY"]),
            allow_redirects=True,
        )
        return response.url


def read_tsv_bytes(f_in: TextIO) -> list[list[str]]:
    """Read TSV from TextIO object and return list of rows (each row a list of str)."""
    reader = csv.reader(f_in, delimiter="\t")
    return [row for row in reader]


def read_tsv(tsv: Path) -> list[list[str]]:
    """Read TSV from Path and return list of rows (each row a list of str)."""
    with tsv.open("rt") as f_in:
        return read_tsv_bytes(f_in)


def sanitize_to_ascii_underscore(text: str) -> str:
    """Replace unicode characters with similar ascii, replace whitespace and commas with underscores."""
    # 1. Normalize Unicode to NFKD form to separate characters from accents
    # 2. Encode to ASCII and ignore characters that cannot be converted
    # 3. Decode back to a string
    text = unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode("ascii")

    # 4. Replace one or more commas or whitespace characters with a single underscore
    # \s+ matches spaces, tabs, and newlines
    return re.sub(r"[,\s]+", "_", text)


def retry[**P, R](
    num_tries: int = 1,
    delay: float = 5.0,
    backoff: float = 2.0,
    no_retry_exceptions: Collection[type] = (),
    logger: logging.Logger | None = None,
) -> Callable[[Callable[P, R]], Callable[P, R]]:
    """Decorator for retrying a stochastic error-prone function.

    Can be invoked directly on undecorated functions like:
    retry(logger=my_logger)(some_func)(func_args_and_kwargs)

    Args:
        num_tries: Total number of attempts (not retries) to make.
        delay: Initial time to wait between retries.
        backoff: Multiplicative factor to increase delay with successive retries.
        no_retry_exceptions: Collection of exceptions that should immediately fail (no retries).
        logger: If supplied, log error messages and retry attempts.
    """
    if num_tries <= 0:
        raise ValueError(
            f"num_tries is total number of attempts, not retries, so must be >= 1. Got: {num_tries}."
        )

    def decorator_repeat(f: Callable[P, R]) -> Callable[P, R]:
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> R:
            nonlocal delay
            nonlocal backoff
            nonlocal logger
            for _ in range(num_tries - 1):
                try:
                    return f(*args, **kwargs)
                except Exception as exception:
                    if type(exception) in no_retry_exceptions:
                        if logger is not None:
                            logger.warning(
                                f"Not retrying exception of type '{type(Exception).__name__}'"
                            )
                        raise
                    if logger is not None:
                        logger.error(f"{exception}")
                        logger.warning(f"Retrying in {delay} seconds...")
                    time.sleep(delay)
                    delay *= backoff
            return f(*args, **kwargs)

        return wrapper

    return decorator_repeat


def download_record(
    record: IgvfRecord,
    igvf_mode: IgvfMode,
    chunk_size: int = 2**20,
    output: Path | None = None,
    logger: logging.Logger | None = None,
) -> Path:
    """Download the file in the provided record, and return the output path.

    Args:
        record: IGVF portal record to download.
        igvf_mode: prod, staging, or sandbox.
        chunk_size: Chunk size for streaming download, in bytes.
        output: If specified, download to that folder using href as file name.
          If unspecified, download in working folder.
        logger: If provided, use for logging.
    Returns:
        Final local path to downloaded file.
    """
    https_url = get_record_http_url(record, igvf_mode=igvf_mode)
    download_path = Path(record["href"].rsplit("/", 1)[-1])
    if output is None:
        output = download_path
    else:
        if (download_path.suffixes) == 0:
            # assume output specified the folder, not the whole path, so fix it
            output = output / download_path.name

    output.parent.mkdir(exist_ok=True, parents=True)

    if logger is not None:
        logger.info(f"Downloading {https_url} to {output}")
    with requests.get(https_url, stream=True) as remote_in:
        remote_in.raise_for_status()
        with output.open("wb") as local_out:
            for chunk in remote_in.iter_content(chunk_size=chunk_size):
                # If you have chunk encoded response uncomment if
                # and set chunk_size parameter to None.
                # if chunk:
                local_out.write(chunk)
    return output


def write_csv[T: FromTypedDict](
    rows: Iterable[T],
    output_csv: Path,
    sep: Literal[",", "\t"] | None = None,
    logger: logging.Logger | None = None,
) -> None:
    """Write rows derived from TypedDict to output CSV."""
    if logger is not None:
        logger.info(f"Writing to {output_csv}")
    if sep is None:
        sep = get_sep_from_path(output_csv)
    output_csv.parent.mkdir(exist_ok=True, parents=True)
    with maybe_gzipped(output_csv, mode="w") as csv_out:
        writer = csv.DictWriter(
            csv_out, fieldnames=T.__annotations__.keys(), delimiter=sep
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def md5sum(filepath: Path, chunk_size: int = 2**20) -> str:
    """Compute md5sum of local file."""
    hasher = hashlib.md5()
    with filepath.open("rb") as f:
        while chunk := f.read(chunk_size):
            hasher.update(chunk)
    return hasher.hexdigest()


# TODO: remove this block once I'm sure it's not needed.
# It requires crcmod as a dependency, so to restore we need to pixi add --pypi crcmod
# # CRC-64/NVME: refin=true, refout=true → rev=True
# # crcmod expects the NORMAL (non-reflected) polynomial when rev=True
# CRC64_NVME_FUNC = crcmod.mkCrcFun(
#     0xAD93D23594C93659,  # normal-form polynomial
#     initCrc=0xFFFFFFFFFFFFFFFF,  # init value
#     xorOut=0xFFFFFFFFFFFFFFFF,  # xorout value (was wrong before!)
#     rev=True,  # refin=true, refout=true
# )
# """Define the CRC64NVME hash function."""

# # TODO: make test
# # # Verify with the known check value: CRC of "123456789" should be 0xAE8B14860A799888
# # check_val = CRC64_NVME_FUNC(b"123456789")
# # assert check_val == 0xAE8B14860A799888, f"Check failed: {hex(check_val)}"
# # print("✅ Check value verified")


# def crc64_nvme_checksum(filepath: Path) -> str:
#     """Compute checksum used by AWS S3."""
#     with filepath.open("rb") as f_in:
#         checksum_int = CRC64_NVME_FUNC(f_in.read())
#     return base64.b64encode(checksum_int.to_bytes(8, "big")).decode()
