import csv
import gzip
import logging
import os
import sys
from collections.abc import (
    Collection,
    Generator,
    Iterator,
)
from contextlib import contextmanager
from io import TextIOWrapper
from pathlib import Path
from typing import TextIO


def setup_logger(logger: logging.Logger, level: int = logging.INFO) -> None:
    logger.handlers.clear()
    logger.setLevel(level)
    ch = logging.StreamHandler(stream=sys.stderr)
    f_formatter = logging.Formatter(
        "%(asctime)s %(name)s %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    ch.setLevel(level)
    ch.setFormatter(f_formatter)
    logger.addHandler(ch)


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


def check_access_keys():
    if "IGVF_API_KEY" not in os.environ or "IGVF_SECRET_KEY" not in os.environ:
        raise ValueError("IGVF_API_KEY and IGVF_SECRET_KEY must be set in environment")


@contextmanager
def maybe_gzipped_reader(maybe_gzipped: Path) -> Generator[TextIO]:
    if maybe_gzipped.name.endswith(".gz"):
        with gzip.open(f"{maybe_gzipped}", "rb") as gzip_in:
            yield TextIOWrapper(gzip_in)
    else:
        yield maybe_gzipped.open("rt")


def iter_csv_rows(
    csv_path, sep: str | None = None, required_columns: Collection[str] | None = None
) -> Iterator[dict[str, str]]:
    if sep is None:
        match csv_path.name.split(".", 1)[-1]:
            case "tsv" | "tsv.gz":
                sep = "\t"
            case "csv" | "csv.gz":
                sep = ","
            case _:
                raise ValueError(
                    f"Unable to automatically determine delimiter for path {csv_path}"
                )
    with maybe_gzipped_reader(csv_path) as f:
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
