from pathlib import Path
from collections.abc import Iterator
import dataclasses

import gzip

from pseudobulk.types import (
    Barcode,
    Contig,
    POS_DTYPE,
)


@dataclasses.dataclass(slots=True, weakref_slot=False, kw_only=True)
class Fragment:
    contig: Contig
    start: POS_DTYPE.type
    end: POS_DTYPE.type
    barcode_sample: Barcode
    num_reads: int

    def __str__(self) -> str:
        return f"{self.contig}\t{self.start}\t{self.end}\t{self.barcode_sample}\t{self.num_reads}\n"

    @classmethod
    def from_line(cls, fragment_line: str) -> "Fragment":
        chro, start, end, barcode_sample, reads = fragment_line.strip().split("\t")
        return cls(
            contig=Contig(chro),
            start=POS_DTYPE.type(start),
            end=POS_DTYPE.type(end),
            barcode_sample=Barcode(barcode_sample),
            num_reads=int(reads),
        )

    @classmethod
    def from_file(cls, fragments_file: Path) -> Iterator["Fragment"]:
        opener = gzip.open if fragments_file.suffix == ".gz" else open
        with opener(f"{fragments_file}", "rt") as fragments_in:
            for line in fragments_in:
                yield cls.from_line(line)

    @property
    def shifted(self) -> "Fragment":
        return dataclasses.replace(self, start=self.start + 4, end=self.end - 4)

    @property
    def start_point(self) -> "Fragment":
        return dataclasses.replace(self, start=self.start, end=self.start + 1)

    @property
    def end_point(self) -> "Fragment":
        return dataclasses.replace(self, start=self.end - 1, end=self.end)
