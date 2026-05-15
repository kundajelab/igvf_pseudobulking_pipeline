# IGVF Pseudobulking

This repository is a tool to assist with IGVF pseudobulking.
It helps you download IGVF single-cell data and perform pseudobulking.

---

## Installation

### 0) Prerequisites

- Ensure you have **download access** to data from the [IGVF portal](data.igvf.org).
- You’ll need an **access key** and **secret key** to the portal.

### 1) Clone this repository

```bash
git clone https://github.com/kundajelab/igvf_pseudobulking_pipeline.git
cd igvf_pseudobulking_pipeline
```

### 2) Install the environment

Use the provided `environment.yml` file to create the `igvf_pseudobulk` conda environment.

```bash
conda env create -f environment.yml
conda activate igvf_pseudobulk
```

### 3) Pick a workspace directory to work in

- A `workspace` is a directory on your filesystem where all downloaded inputs, intermediate files, and outputs will be placed. You choose this path; the pipelines will create and organize subfolders inside it.

Set your workspace path (example):

```bash
workspace=/path/to/dir
```

---

## Usage

You can run the IGVF pseudobulking code via:

```bash
python igvf_process.py <download|pseudobulk|full> [options]
```

### Download pipeline

The `download` option will download the uniformly processed fragments and RNA counts for the specified IGVF analysis accessions. You can run it with:

```bash
python igvf_process.py download \
 -al ${analysis_accession_list} \
 -w ${workspace} \
 -ak ${access_key} \
 -sk ${secret_key}
```

The arguments are:

- `-al` - a comma-separated list of IGVF analysis accessions to download data from (ex: IGVFDS6430MYNQ,IGVFDS5417HJRJ)
- `-w` - your chosen [workspace directory](#3-pick-a-workspace-directory-to-work-in)
- `ak` - your access key
- `sk` - your secret key

Once it is done running, you will find the following directories inside your `workspace`:

- `raw_fragments/` — this directory will contain one fragments file named `${ANALYSIS_ACCESSION}.tsv.gz` per IGVF analysis accession will be inside it
- `raw_rna/` - this directory will contain one AnnData RNA count matrix named `${ANALYSIS_ACCESSION}.h5ad` per IGVF analysis accession will be inside it

### Pseudobulk pipeline

The `pseudobulk` option uses a standardized cell annotation file to aggregate cells into pseudobulks and perform QC. You can run it with:

```bash
python igvf_process.py pseudobulk \
 -w ${workspace} \
 -a /path/to/annotation.tsv \
 -s ${species} \
 -c ${numCPUs}
```

The required arguments are:

- `-w` - your chosen [workspace directory](#3-pick-a-workspace-directory-to-work-in)
- `-a` - path to your annotation file (see: [Annotation File Requirements](#annotation-file-requirements))
- `-s` - the species (must be one of `human` or `mouse`) reference files are sourced automatically from `genome_data/${species}`
- `-c` - number of CPUs to use for parallelization

Once it is done running, you will find the following inside your `workspace`:

1) `cell_name_to_annotation_mapping.tsv` — This file contains a mapping from `cell_name` values in the annotation file to a unique `annotation_ID` of the form `annotation_NUMBER`. This internal mapping is done so no restrictions are placed on the format of `cell_name`. You can refer to this file in order to identify which pseudobulks correspond to which cell names.
2) `pseudobulk/` — This directory will contain one subdirectory per pseudobulk. These will be named `${annotation_ID}-${subsample}` where `annotation_ID` is the unique identifier of a cell name. This can be deciphered by referencing `cell_name_to_annotation_mapping.tsv`. Each pseudobulk subdirectory contains 7 files:
    - `fragments.tsv` — combined fragments from all cells in the pseudobulk.
    - `raw_insertions.bw` — pileup of +4/-4 shifted insertions from all fragments.
    - `peaks.narrowPeak` — peaks called on shifted insertions (reproducible peak calling with MACS3).
    - `peaks_minuslog10pval.bw` — -log10 p-value signal track produced by MACS3.
    - `rna_counts_mtx.h5ad` — RNA counts matrix subset to cells in the pseudobulk.
    - `pseudobulk_expression.tsv` — pseudobulked expression values.
    - `per_cell_qc.tsv` — per-cell QC for cells in the pseudobulk with columns:
      - `analysis_accession` - the analysis accession the cell came from
      - `barcode` - the cell's barcode
      - `subsample` - the subsample the cell derives from
      - `rna_read_count` - the number of RNA reads
      - `gene_count` - the number of unique genes
      - `pct_mito` - the percent of reads from mitochondrial genes
      - `pct_ribo` - the percent of reads from ribosomal genes
      - `num_frags` - the number of ATAC fragments
      - `pct_duplicated_reads` - the percent of ATAC reads that are duplicates
      - `nucleosomal_signal` - the fraction of mono-nucleosomal fragments (148-294 bp) to the fraction of nucleosome-free fragments (1-147 bp)
      - `tss_enrichment` - the enrichment of shifted insertions near TSS centers over the average number of shifted insertions at TSS flanks
      - `frip` - the fraction of fragments overlapping peaks (computed on the pseudobulk peaks)
3) `analysis_accession_qc_reports/` — this directory will contain one QC file per analysis accession named `${analysis_accession}_per_cell_qc.tsv`. This QC file contains one row per barcode that appeared in the analysis accession fragment file or RNA counts matrix. Not every barcode is guaranteed to be a real annotated cell, however. The QC columns are:
    - `analysis_accession` - the analysis accession the barcode came from
    - `barcode` - the barcode
    - `annotated` - whether or not this barcode was annotated; may NOT be true for all barcode in this file
    - `found_in_rna` - whether or not this barcode was found in the RNA counts matrix
    - `found_in_atac` - whether or not this barcode was found in the ATAC fragments file
    - `rna_read_count` - the number of RNA reads
    - `gene_count` - the number of unique genes
    - `pct_mito` - the percent of reads from mitochrondrial genes
    - `pct_ribo` - the percent of reads from ribosomal genes
    - `num_frags` - the number of ATAC fragments
    - `pct_duplicated_reads` - the percent of ATAC reads that are duplicates
    - `nucleosomal_signal` - the fraction of mono-nucleosomal fragments (148-294 bp) to the fraction of nucleosome-free fragments (1-147 bp)
    - `tss_enrichment` - the enrichment of shifted insertions near TSS centers over the average number of shifted insertions at TSS flanks
4) `pseudobulk_qc.tsv` - this file contains per-pseudobulk QC metrics. The QC metrics are not the mean or median of the corresponding QC metrics from the individual cells in the pseudobulk. Rather, these QC metrics are computed when treating the cluster as a single pseudobulked entity. The columns in this file are:
    - `directory_name` - the name of the subdirectory in the `pseudobulk/` directory that corresponds to this pseudobulk in the form of `{annotation_ID}-{subsample}`
    - `pseudobulk` - the name of the pseudobulk in the form of `{cell_name}-{subsample}`
    - `cell_name` - the cell name itself
    - `subsample` - the subsample
    - `num_cells` - the number of cells in the pseudobulk
    - `rna_read_count` - the number of RNA reads
    - `gene_count` - the number of unique genes
    - `pct_mito` - the percent of reads from mitochrondrial genes
    - `pct_ribo` - the percent of reads from ribosomal genes
    - `num_frags` - the number of ATAC fragments
    - `pct_duplicated_reads` - the percent of ATAC reads that are duplicates
    - `nucleosomal_signal` - the fraction of mono-nucleosomal fragments (148-294 bp) to the fraction of nucleosome-free fragments (1-147 bp)
    - `tss_enrichment` - the enrichment of shifted insertions near TSS centers over the average number of shifted insertions at TSS flanks
    - `frip` - the fraction of fragments overlapping peaks (computed on the pseudobulk peaks)

### Full pipeline

The `full` option runs the Download pipeline followed by the Pseudobulk pipeline, inferring which analysis accessions to download from the provided annotation file. You can run it as

```bash
python igvf_process.py full \
 -w ${workspace} \
 -a /path/to/annotation.tsv \
 -s ${species} \
 -ak ${access_key} \
 -sk ${secret_key} \
 -c ${numCPUs}
```

The required arguments are:

- `-w` - your chosen [workspace directory](#3-pick-a-workspace-directory-to-work-in)
- `-a` - path to your annotation file (see: [Annotation File Requirements](#annotation-file-requirements))
- `-s` - species, one of `human` or `mouse`; reference files are sourced automatically from `genome_data/<species>`
- `-ak` - your access key
- `-sk` - your secret key
- `-c` - number of CPUs to use for parallelization

The outputs are the outputs of the Download and Pseudobulk pipelines as described above.
Outputs are the union of the Download and Pseudobulk results described above, created within your `workspace`.

---

## Annotation File Requirements

Annotations are expected to be tab-separated files with the following properties:

- Required columns:
  - `barcode_sample`: The barcode, exactly as you would find it in the fragments file.
  - `cell_name`: **The primary annotation upon which pseudobulking is performed**.
  - `cell_description`: A free text description unique to the cell name (does not affect pseudobulking).
  - `CL_id`: The CLID closest to the cell name (does not affect pseudobulking).
  - `CL_term_name`: The term name of the CL_id (does not affect pseudobulking).
  - `subsample`: The subsample from which the cell came from. **Pseudobulking occurs at the `cell_name` x `subsample` level**.
  - `analysis_set_accession`: The accession of the analysis set. This cell's barcode_sample will be looked for in the fragment file corresponding to the `analysis_set_accession`.
- Annotation file format requirements:
  - The annotation file **must be a .tsv file**.
  - No column values should contain tabs (otherwise, file parsing will fail).
  - The `subsample` column cannot contain hyphens.

---

## Genome Data

The `genome_data/` directory holds reference data for supported species (`human`, `mouse`). It includes chromosome sizes, blacklist regions, gene information, and TSS data used during processing.

For details on how these files were generated, see `genome_data/README.md`.

## Processing and calculation notes

- Fragments that do not come from a chromosome in the speciesl chromosome sizes file are **ignored entirely**. They will not even contribute towards the `num_frags` QC metric.
- The `tss_enrichment` is computed per barcode via the following procedure:
  - For each barcode, a generic window of 2001 bases (center +- 1000 bases) around TSS centers is instantiated. Index 0 corresponds to -1000 bases from all TSSs and index 2000 corresponds to +1000 bases from all TSSs.
  - From each ATAC fragment, two shifted (+4/-4) insertions are generated.
  - For each insertion, its distance to the closest TSS (based on the species TSS file) is computed. (chrM TSSs are ignored)
  - If the distance is within 1000 bp, the signed distance is computed (strand is taken into account). The value in the barcode's generic TSS window corresponding to that signed distance is incremented, representing finding an insertion at that distance from a TSS.
  - After all fragments have been processed, the average number of insertions per base in the flanking 200 bases of the generic TSS window are computed as `tss_insertions_flank_mean = (np.sum(tss_insertions[:100]) + np.sum(tss_insertions[-100:]))/200`.
  - Then, the average number of insertions per base in the central 5 bases are computed as `tss_insertions_center = np.mean(tss_insertions[1000-2:1000+3])`.
  - The TSS Enrichment is computed as the ratio between these values: `(tss_insertions_center/(tss_insertions_flank_mean if tss_insertions_flank_mean > 0 else 1))`.
- When computing `nucleosomal_signal`, the fragment size is computed based on the +4/-4 shifted insertion fragment endpoints rather than the raw fragment endpoints. Also, a pseudocount of 1 is applied to both mono-nucleosomal fragments (148-294 bp) and nucleosome-free fragments (1-147 bp).
