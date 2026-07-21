from collections.abc import (
    Collection,
    Iterator,
)
import csv
from io import TextIOWrapper
from typing import (
    Generator,
    TextIO,
)
from pathlib import Path
from contextlib import contextmanager
import gzip


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
        for row in reader:
            yield row


def iter_pseudobulk_dirs(pseudobulk_dir: Path) -> Iterator[Path]:
    """List all the sub-folders in the pseudobulk_dir. Each should correspond to a unique pseudobulk."""
    for folder in sorted(pseudobulk_dir.iterdir()):
        if folder.is_dir():
            yield folder
