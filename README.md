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

The `pseudobulk` option aggregates cells into pseudobulks based on annotation columns and generates per-pseudobulk data and QC. You can run it with:

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

You can also use the optional flag:

- `--at_annotation_level` - whether or not to pseudobulk at the annotation level. If not set (the default behavior), pseudobulks will be generated at the annotation x subsample level. If set, pseudobulks will be generated at the annotation level and pooled across subsamples.

Once it is done running, you will find the following inside your `workspace`:

1) `pseudobulk/` — this directory will contain one subdirectory per pseudobulk. By default, these will be named `${annotation-column}-${annotation-value}-${subsample}`. If the `--at_annotation_level` flag is used, the pseudobulks will be named `${annotation-column}-${annotation-value}`.Each pseudobulk subdirectory contains 7 files:
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
2) `analysis_accession_qc_reports/` — this directory will contain one QC file per analysis accession named `${analysis_accession}_per_cell_qc.tsv`. This QC file contains one row per barcode that appeared in the analysis accession fragment file or RNA counts matrix. Not every barcode is guaranteed to be a real annotated cell, however. The QC columns are:
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
3) `pseudobulk_qc.tsv` - this file contains per-pseudobulk QC metrics. The QC metrics are not the mean or median of the corresponding QC metrics from the individual cells in the pseudobulk. Rather, these QC metrics are computed when treating the cluster as a single pseudobulked entity. The columns in this file are:
    - `pseudobulk` - the name of the pseudobulk
    - `annotation_level` - the annotation level (the column in the annotation file from which this pseudobulk was derived from)
    - `annotation` - the annotation itself
    - `subsample` - the subsample. **NOTE: THIS COLUMN WILL NOT BE PRESENT IF THE `--at_annotation_level` flag is used.**
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

You can also use the optional flag:

- `--at_annotation_level` - whether or not to pseudobulk at the annotation level. If not set (the default behavior), pseudobulks will be generated at the annotation x subsample level. If set, pseudobulks will be generated at the annotation level and pooled across subsamples.

The outputs are the outputs of the Download and Pseudobulk pipelines as described above.
Outputs are the union of the Download and Pseudobulk results described above, created within your `workspace`.

---

## Annotation File Requirements

Annotations are expected to be tab-separated files with the following properties:

- Required columns: `analysis_accession`, `barcode`, `subsample`, `annotation`.
- `annotation` is considered to be the primary cell annotation. However, multiple levels of annotations can be specified (which may be relevant for cell subtypes). Any other column beginning with `annotation` will be considered to be another level of annotation and pseudobulks will be made for it.
- the `subsample` column and no `annotation` column can contain spaces, hyphens (`-`), or other nonstandard characters.

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
