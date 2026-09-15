import logging
import re
from collections.abc import Iterable, Iterator
from pathlib import Path

from igvf_portal import utils
from igvf_portal.connection import PConnection
from igvf_portal.constants import VERSION
from igvf_portal.enums import IgvfMode
from igvf_portal.types import (
    AccessionId,
    GeneInfoRow,
    IgvfRecord,
    PortalId,
    TssRow,
)

_ID_MATCH = re.compile(r'gene_id "([^"]+)"')
_NAME_MATCH = re.compile(r'gene_name "([^"]+)"')


def _get_reference_file_ids(
    connection: PConnection, file_accessions: Iterable[AccessionId]
) -> set[PortalId]:
    return {
        ref if isinstance(ref, str) else ref["@id"]
        for file_accession in file_accessions
        for ref in connection.lookup_record(file_accession, frame="object")["reference_files"]
    }


def _get_reference_file_records(
    connection: PConnection, reference_file_ids: Iterable[PortalId]
) -> tuple[IgvfRecord, IgvfRecord]:
    # content_type, href and s3_uri are all the file's own properties
    reference_file_records = [
        connection.lookup_record(reference_file_accession, frame="object")
        for reference_file_accession in reference_file_ids
    ]
    match [_rec for _rec in reference_file_records if _rec["content_type"] == "genome reference"]:
        case [fasta_reference]:
            pass
        case _ as fasta_references:
            raise ValueError(f"Expected 1 fasta reference, got: {len(fasta_references)}")
    match [
        _rec for _rec in reference_file_records if _rec["content_type"] == "transcriptome reference"
    ]:
        case [transcriptome_reference]:
            pass
        case _ as transcriptome_references:
            raise ValueError(
                f"Expected 1 transcriptome reference, got: {len(transcriptome_references)}"
            )

    return (fasta_reference, transcriptome_reference)


def _iter_gtf_rows(gtf_path: Path, logger: logging.Logger) -> Iterator[list[str]]:
    """Yield lines of the gene_map to output."""
    logger.info(f"Parsing GTF: {gtf_path}")
    with utils.maybe_gzipped(gtf_path, mode="r") as gtf_in:
        for line in gtf_in:
            if line.startswith("#"):
                continue
            yield line.strip().split("\t")


def _iter_gene_map(gtf_path: Path, logger: logging.Logger) -> Iterator[GeneInfoRow]:
    """Yield lines of the gene_map to output."""
    for gtf_row_parts in _iter_gtf_rows(gtf_path=gtf_path, logger=logger):
        if gtf_row_parts[2] != "gene":
            continue

        attributes = gtf_row_parts[8]
        id_match = _ID_MATCH.search(attributes)
        name_match = _NAME_MATCH.search(attributes)

        if id_match is not None:
            gene_id = id_match.group(1)
            # Strip version if your matrix doesn't have them
            # gene_id = gene_id.split('.')[0]

            gene_name = name_match.group(1) if name_match else gene_id

            # Check Mito/Ribo status HERE
            is_mito = gene_name.startswith(("MT-", "mt-"))
            is_ribo = gene_name.startswith(("RPS", "RPL", "Rps", "Rpl"))
            yield {
                "gene_id": gene_id,
                "gene_name": gene_name,
                "mt": is_mito,
                "ribo": is_ribo,
            }


def _iter_tss_rows(gtf_path: Path, logger: logging.Logger) -> Iterator[TssRow]:
    """Yield Transcription Start Site (TSS) rows for transcripts in GTF."""
    for line_split in _iter_gtf_rows(gtf_path=gtf_path, logger=logger):
        chro, feature, start, end, strand, info = (
            line_split[0],
            line_split[2],
            int(line_split[3]),
            int(line_split[4]),
            line_split[6],
            line_split[8],
        )
        if feature == "transcript":
            info_dict = {
                key: val.strip('"') for key, val in (x.split(" ", 1) for x in info.split("; "))
            }
            if strand not in ("+", "-"):
                raise ValueError(f"Invalid strand: {strand}")
            yield {
                "gene": info_dict["gene_name"],
                "transcript": info_dict["transcript_name"],
                "chro": chro,
                "TSS": start if strand == "+" else end,
                "strand": strand,
            }


def get_references(
    key: str,
    *,
    igvf_mode: IgvfMode = IgvfMode.prod,
    output: Path | None = None,
    chunk_size: int = 2**20,
    gene_info_name: str = "gene_info.csv",
    tss_name: str = "tss.tsv",
) -> None:
    """Download the reference FASTA and GTF used by file(s) described by key.

    Then form gene_info.csv and tss.tsv from the GTF.

    Args:
        key: Comma-separated list of alias or accession ID to look for reference file
            dependencies.
        igvf_mode: Mode for accessing the IGVF Portal.
        output: If specified, download to that folder and write gene info and TSS files there. If
            unspecified, use the working folder.
        chunk_size: Chunk size for streaming download, in bytes.
        gene_info_name: File name of the gene info CSV to write.
        tss_name: File name of the TSS TSV to write.
    """
    utils.check_access_keys()
    logger = utils.get_logger_from_file(__file__)
    logger.info(f"Version: {VERSION}")

    connection = PConnection.new(igvf_mode=igvf_mode)
    split_keys = {AccessionId(_split_key.strip()) for _split_key in key.split(",")}
    reference_file_accessions = _get_reference_file_ids(
        connection=connection, file_accessions=split_keys
    )
    fasta_record, gtf_record = _get_reference_file_records(
        connection=connection, reference_file_ids=reference_file_accessions
    )

    _ = connection.download_record(fasta_record, chunk_size=chunk_size, output=output)
    gtf_path = connection.download_record(gtf_record, chunk_size=chunk_size, output=output)

    if output is None:
        output = Path(".")

    # write gene_info CSV
    gene_map_rows_iter = _iter_gene_map(gtf_path=gtf_path, logger=logger)
    utils.write_csv(
        rows=gene_map_rows_iter,
        output_csv=output / gene_info_name,
        row_type=GeneInfoRow,
        logger=logger,
    )
    # write tss TSV
    tss_rows_iter = _iter_tss_rows(gtf_path=gtf_path, logger=logger)
    utils.write_csv(
        rows=tss_rows_iter, output_csv=output / tss_name, row_type=TssRow, logger=logger
    )
