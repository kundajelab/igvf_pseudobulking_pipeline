import dataclasses
import logging
from collections.abc import Iterable
from enum import StrEnum
from pathlib import Path
from typing import Final, NewType

import numpy as np
import scipy.sparse

Barcode = NewType("Barcode", str)
Contig = NewType("Contig", str)
PseudobulkName = NewType("PseudobulkName", str)

POS_DTYPE: Final[np.dtype] = np.dtype(np.int32)
type POS_ARRAY = np.ndarray[tuple[int], np.dtype[POS_DTYPE.type]]
COUNTS_DTYPE: Final[np.dtype] = np.dtype(np.uint16)
type COUNTS_ARRAY = np.ndarray[tuple[int], np.dtype[COUNTS_DTYPE.type]]
type COUNTS_MATRIX = scipy.sparse.csr_array[COUNTS_DTYPE.type, tuple[int, int]]

_REPO_DIR: Final[Path] = Path(__file__).parent.parent.parent.parent
_GENOME_DATA_DIR: Final[Path] = _REPO_DIR / "genome_data"


class FailureAction(StrEnum):
    """Enum for what action to take when a check for expected values fails."""

    exception = "exception"
    """raise an exception."""

    warning = "warning"
    """log a warning."""

    sentinal = "sentinal"
    """create a sentinal file."""


@dataclasses.dataclass
class FailureHandler:
    """Takes the configured actions when a check for expected values fails."""

    actions: Iterable[FailureAction]
    sentinal_file_name: Path
    logger: logging.Logger

    def handle_failure(self, message: str) -> None:
        """Take each of the configured actions for a failed check.

        Args:
            message: Description of the failed check.
        """
        for action in self.actions:
            match action:
                case FailureAction.exception:
                    raise RuntimeError(message)
                case FailureAction.warning:
                    self.logger.warning(message)
                case FailureAction.sentinal:
                    self.sentinal_file_name.parent.mkdir(parents=True, exist_ok=True)
                    with open(self.sentinal_file_name, mode="a") as f_out:
                        f_out.write(f"{message}\n")
