import csv
import dataclasses
import gzip
import hashlib
import logging
import os
import re
import sys
import time
import unicodedata
from collections.abc import (
    Callable,
    Collection,
    Generator,
    Iterable,
    Iterator,
    Mapping,
)
from contextlib import contextmanager
from io import TextIOWrapper
from pathlib import Path
from typing import Literal, TextIO

import requests

from igvf_portal.types import IgvfRecord


def setup_logger(logger: logging.Logger, level: int = logging.INFO) -> None:
    """Remove all handlers from `logger` and set its level."""
    logger.handlers.clear()
    logger.setLevel(level)


def get_logger(name: str | None = None, level: int = logging.INFO) -> logging.Logger:
    """Get the named logger (root logger if name is None), cleared of handlers and set to level."""
    logger = logging.getLogger(name)
    setup_logger(logger, level)
    return logger


def get_logger_from_file(file_name: str, level: int = logging.INFO) -> logging.Logger:
    """Get a logger named "<package> <tool>" from the path of a tool module (i.e. `__file__`)."""
    # file_name is <...>/<package>/tools/<tool>.py. Splitting the string on "." would break on any
    # dotted directory above it, e.g. lib/python3.14/site-packages
    tool_path = Path(file_name)
    package_name = tool_path.parent.parent.name
    tool_name = tool_path.stem
    return get_logger(
        f"{package_name.replace('_', '-')} {tool_name.replace('_', '-')}", level=level
    )


def fix_igvf_logging(level: int = logging.INFO, debug_logger: logging.Logger | None = None) -> None:
    """Configure logging for IGVF Portal tools.

    Sends root logger output to stderr, silences the "iu_error" logger, and keeps the AWS SDK
    loggers no more verbose than INFO.

    Args:
        level: Log level for the root logger and debug logger.
        debug_logger: Logger to set up at `level`. Defaults to the "iu_debug" logger.
    """
    root = get_logger(name=None, level=level)
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
    # At DEBUG the AWS SDK logs every signed S3 request with its headers, including the temporary
    # session token, and a large upload makes tens of thousands of requests. Keep it at INFO, or at
    # the requested level if that is quieter.
    for aws_logger_name in ("botocore", "boto3", "s3transfer"):
        logging.getLogger(aws_logger_name).setLevel(max(level, logging.INFO))


def check_access_keys() -> None:
    """Raise ValueError unless IGVF Portal API keys are set in the environment."""
    if "IGVF_API_KEY" not in os.environ or "IGVF_SECRET_KEY" not in os.environ:
        raise ValueError("IGVF_API_KEY and IGVF_SECRET_KEY must be set in environment")


_HTTP_TIMEOUT: float = 60.0
"""Timeout in seconds for HTTP requests made outside of a PConnection."""


def lookup_ontology_by_cl_id(cl_id: str) -> list[dict[str, object]] | None:
    """Retrieve the full OLS record for a term given its CL_id or short_form ID.

    Example: CL:0011026
    """
    response = requests.get(
        "https://www.ebi.ac.uk/ols4/api/terms",
        params={"short_form": cl_id.replace(":", "_")},
        timeout=_HTTP_TIMEOUT,
    )
    response.raise_for_status()
    data = response.json()

    # The response is a paginated object; terms are in the embedded "terms" list
    return data.get("_embedded", {}).get("terms", None)


@contextmanager
def maybe_gzipped(maybe_gzipped: Path, mode: Literal["w", "r"]) -> Generator[TextIO]:
    """Get a text reader from a path that may or may not be gzipped.

    When writing, set mtime to 0 to ensure deterministic output.
    """
    # NOTE: the TextIOWrapper must be closed (flushing its buffer) before the GzipFile it wraps,
    # otherwise the buffered text never reaches the gzip stream and the output is truncated.
    if maybe_gzipped.name.endswith(".gz"):
        if mode == "r":
            with (
                gzip.open(f"{maybe_gzipped}", "rb") as f_gzip,
                TextIOWrapper(f_gzip) as f_text,
            ):
                yield f_text
        else:
            with (
                maybe_gzipped.open("wb") as raw_in,
                gzip.GzipFile(mode="wb", fileobj=raw_in, mtime=0) as f_gzip,
                TextIOWrapper(f_gzip) as f_text,
            ):
                yield f_text
    else:
        with maybe_gzipped.open("rt" if mode == "r" else "wt") as f_text:
            yield f_text


def get_sep_from_path(csv_or_tsv: Path) -> Literal["\t", ","]:
    """Guess the file separator from the path."""
    match csv_or_tsv.suffixes:
        case [*_parts, ".tsv"] | [*_parts, ".tsv", ".gz"]:
            return "\t"
        case [*_parts, ".csv"] | [*_parts, ".csv", ".gz"]:
            return ","
        case _:
            raise ValueError(
                f"Could not infer separator for file {csv_or_tsv} with suffix "
                f"{csv_or_tsv.suffix}. Please specify separator with 'sep' argument."
            )


def iter_csv_rows(
    csv_path: Path, sep: str | None = None, required_columns: Collection[str] | None = None
) -> Iterator[dict[str, str]]:
    """Iterate over rows of a (possibly gzipped) CSV or TSV file as dicts.

    Args:
        csv_path: Path to the CSV or TSV file.
        sep: Column separator. If None, infer it from the file suffix.
        required_columns: If supplied, raise ValueError if any of these columns are missing.
            Whitespace is stripped from column names before checking.

    Yields:
        Each row as a dict mapping column name to value.
    """
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
                raise ValueError(f"{csv_path} is missing required columns: {','.join(missing)}")
        yield from reader


def iter_pseudobulk_dirs(pseudobulk_dir: Path) -> Iterator[Path]:
    """List all the sub-folders in the pseudobulk_dir.

    Each should correspond to a unique pseudobulk.
    """
    for folder in sorted(pseudobulk_dir.iterdir()):
        if folder.is_dir():
            yield folder


def parse_s3_uri(s3_uri: str) -> tuple[str, str]:
    """Get bucket and object name from AWS S3 URI."""
    if not s3_uri.startswith("s3://"):
        raise ValueError(f"'{s3_uri}' is not a valid S3 path.")
    _, _, bucket_name, object_key = s3_uri.split("/", 3)
    return bucket_name, object_key


def read_tsv_bytes(f_in: TextIO) -> list[list[str]]:
    """Read TSV from TextIO object and return list of rows (each row a list of str)."""
    reader = csv.reader(f_in, delimiter="\t")
    return list(reader)


def read_tsv(tsv: Path) -> list[list[str]]:
    """Read TSV from Path and return list of rows (each row a list of str)."""
    with tsv.open("rt") as f_in:
        return read_tsv_bytes(f_in)


def sanitize_to_ascii_underscore(text: str) -> str:
    """Replace unicode with similar ascii.

    Also replace periods, commas, slashes, and whitespace with underscores.
    """
    # 1. Normalize Unicode to NFKD form to separate characters from accents
    # 2. Encode to ASCII and ignore characters that cannot be converted
    # 3. Decode back to a string
    text = unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode("ascii")

    # 4. Replace one or more periods, commas, forward or back slashes, or whitespace characters with
    #    a single underscore
    # \s+ matches spaces, tabs, and newlines
    return re.sub(r"[-.,/\\\s]+", "_", text)


@dataclasses.dataclass(frozen=True, kw_only=True, slots=True)
class RetryPolicy:
    """How often, and how patiently, to retry a failing operation (see retry)."""

    num_tries: int = 3
    """Total number of attempts (not retries) to make."""
    delay: float = 5.0
    """Initial time in seconds to wait between attempts."""
    backoff: float = 2.0
    """Multiplicative factor to increase delay with successive attempts."""


def retry[**P, R](
    num_tries: int = 1,
    delay: float = 5.0,
    backoff: float = 2.0,
    no_retry_exceptions: Collection[type] = (),
    should_retry: Callable[[Exception], bool] | None = None,
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
        should_retry: If supplied, an exception for which it returns False immediately fails
            (no retries), e.g. to tell a permanent HTTP error from a temporary one.
        logger: If supplied, log error messages and retry attempts.
    """
    if num_tries <= 0:
        raise ValueError(
            "num_tries is total number of attempts, not retries, so must be >= 1. "
            f"Got: {num_tries}."
        )

    no_retry_types = tuple(no_retry_exceptions)

    def decorator_repeat(f: Callable[P, R]) -> Callable[P, R]:
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> R:
            # each call starts from the initial delay, rather than where the last call left off
            call_delay = delay
            for _ in range(num_tries - 1):
                try:
                    return f(*args, **kwargs)
                except Exception as exception:
                    if isinstance(exception, no_retry_types) or (
                        should_retry is not None and not should_retry(exception)
                    ):
                        if logger is not None:
                            logger.warning(
                                f"Not retrying exception of type '{type(exception).__name__}'"
                            )
                        raise
                    if logger is not None:
                        logger.error(f"{exception}")
                        logger.warning(f"Retrying in {call_delay} seconds...")
                    time.sleep(call_delay)
                    call_delay *= backoff
            return f(*args, **kwargs)

        return wrapper

    return decorator_repeat


def write_csv[T: Mapping[str, object]](
    rows: Iterable[T],
    output_csv: Path,
    row_type: type[T],
    sep: Literal[",", "\t"] | None = None,
    logger: logging.Logger | None = None,
) -> None:
    """Write rows derived from TypedDict row_type to output CSV.

    row_type supplies the column names: the type parameter T does not exist at runtime.
    """
    if logger is not None:
        logger.info(f"Writing to {output_csv}")
    if sep is None:
        sep = get_sep_from_path(output_csv)
    output_csv.parent.mkdir(exist_ok=True, parents=True)
    with maybe_gzipped(output_csv, mode="w") as csv_out:
        writer = csv.DictWriter(csv_out, fieldnames=row_type.__annotations__.keys(), delimiter=sep)
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


def file_not_in_portal(record: IgvfRecord) -> bool:
    """Whether the file in `record` was deleted or its upload was not found on the Portal."""
    return record["status"] == "deleted" or record.get("upload_status", "") == "file not found"
