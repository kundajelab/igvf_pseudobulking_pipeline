import logging
import warnings
from collections.abc import Callable

# google-crc32c has no prebuilt C extension for the free-threaded CPython
# build we use, so it always falls back to its pure-Python implementation.
# Filter before importing the tools below, since some import google_crc32c
# transitively and it warns at import time.
warnings.filterwarnings(
    "ignore", message="As the c extension couldn't be imported", category=RuntimeWarning
)

import defopt  # noqa: E402

from pseudobulk.tools.combine_accession_qc import combine_accession_qc  # noqa: E402
from pseudobulk.tools.pseudobulk_rna import pseudobulk_rna  # noqa: E402
from pseudobulk.tools.split_fragments import split_fragments  # noqa: E402
from pseudobulk.tools.summarize_pseudobulk_qc import summarize_pseudobulk_qc  # noqa: E402

tools: list[Callable] = [
    split_fragments,
    pseudobulk_rna,
    summarize_pseudobulk_qc,
    combine_accession_qc,
]


def main() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    defopt.run(tools)


if __name__ == "__main__":
    main()
