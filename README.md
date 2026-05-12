# IGVF Pseudobulking

This repository is a tool to assist with IGVF pseudobulking.
It downloads IGVF single-cell data and performs pseudobulking.

---

## Contents:
- [Installation](#installation)
- [Usage](#usage)
- [Outputs](#outputs)
- [Annotation File Requirements](#annotation-file-requirements)
- [Genome Data](#genome-data)
- [Processing and calculation notes](#processing-and-calculation-notes)
- [Developing](#developing)


## Installation

### 0) Prerequisites

- Ensure you have **download access** to data from the [IGVF portal](data.igvf.org).
- You’ll need an **access key** and **secret key** to the portal.
- Software pre-reqs:
  - git (to clone)
  - pixi (to install and run the pipeline): https://pixi.prefix.dev/latest/installation/
  - the Pseudobulk environment requires
    - gcc
    - openblas
    - hdf5
  - gcc and openblas (on sherlock: `module load system curl/8.17.0 ncurses/6.4 gcc/12.4.0 openblas/0.3.28 hdf5/1.14.4`)

### 1) Clone this repository

```bash
git clone https://github.com/kundajelab/igvf_pseudobulking_pipeline.git
cd igvf_pseudobulking_pipeline
```

### 2) Install the environment

The main project environment is managed by pixi. There are also two sub-projects: `pseudobulk`
(managed by `uv`) and `download_accession_files` (managed by `pixi`). You can just run commands, but
to install up front, open a terminal in the project base folder and run

```bash
pixi run install-all
```

On an HPC cluster like `sherlock`, you may need to launch a compute node (e.g. via `sh_dev`) to
install. On the login nodes, sometimes the install fails due to not enough available file
descriptors.

In the following commands, you can avoid the need for starting commands with  `pixi run` by
activating the environment:
```bash
eval $(pixi shell-hook)
```
This README will keep the pixi run invocations to make clear when pixi is launching commands.

Save your secrets securely to be used by nextflow
```bash
pixi run nextflow secrets set IGVF_API_KEY [YOUR_ACCESS_KEY]
pixi run nextflow secrets set IGVF_SECRET_KEY [YOUR_SECRET_KEY]
```

## Usage

More most use-cases, you can run the IGVF pseudobulking pipeline by opening a terminal in the
igvf_pseudobulking_pipeline project folder and running:

```bash
pixi run pipeline PATH_TO_PSEUDOBULKING_METDATA_TSV
```

This will rapidly download all needed files in parallel and run the pipeline locally. You can also
override the defaults for queue, profile, and output folder if needed (see the `pipeline` task in
[pixi.toml](./pixi.toml)). e.g. to run in the owners queue and change the output folder:
```bash
pixi run pipeline PATH_TO_PSEUDOBULKING_METDATA_TSV owners slurm,apptainer "$HOME/my/output/folder"
```
Or if even more customizability is needed, you can launch nextflow directly without a helper task:
```bash
pixi run nextflow . --metadata_file [OTHER_ARGS]
```
To see the complete list of overrideable parameters, see [nextflow.config](./nextflow.config) (_params_ section).
Although this way you will at a minimum need to manually specify the profile and output folders.

You can find a report with details of execution duration, resource usage, etc from the latest run by
running
```bash
pixi run get-latest-report
```
It will display an HTML file that can be opened in any web browser.

Once you are done, you can delete all the output and temporary files created by the pipeline by
running
```bash
pixi run clear
```

## Outputs
Once the pipeline is done running, you will find the following directories inside your `workspace`:

- `raw_fragments/` — this directory will contain one fragments file named
  `${ANALYSIS_ACCESSION}.tsv.gz` per IGVF analysis accession will be inside it
- `raw_rna/` - this directory will contain one AnnData RNA count matrix named
  `${ANALYSIS_ACCESSION}.h5ad` per IGVF analysis accession will be inside it
- `trace/` - this directory contains information about the pipeline resource usage and execution
  sequence for past runs.
- `work/` - this directory contains nested subfolders with intermediate files of nextflow processes.
  It should not be removed until you are done with outputs, as outputs are generally symlinked from
  within this folder.
- `output/` - this directory contains final output files of the pseudobulking pipeline.
  -  `pseudobulks/` — this directory will contain one subdirectory per pseudobulk. These will be
     named `annotation_${annotation-index}-${subsample}`. Each pseudobulk subdirectory contains 7
     files:
      - `fragments.tsv.gz` — combined fragments from all cells in the pseudobulk.
      - `raw_insertions.bw` — pileup of +4/-4 shifted insertions from all fragments.
      - `peaks.narrowPeak.gz` — TSV with peaks called on shifted insertions (reproducible peak calling with MACS3).
      - `peaks_minuslog10pval.bw` — -log10 p-value signal track produced by MACS3.
      - `rna_counts_mtx.h5ad` — RNA counts matrix subset to cells in the pseudobulk.
      - `pseudobulk_expression.tsv.gz` — pseudobulked expression values.
      - `per_cell_qc.tsv.gz` — per-cell QC for cells in the pseudobulk with columns:
        - `analysis_set_accession` - the analysis accession the cell came from
        - `barcode_sample` - the cell's barcode
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
  - `analysis_accession_qc_reports/` — this directory will contain one QC file per analysis accession
    named `${analysis_accession}_per_cell_qc.tsv.gz`. This QC file contains one row per barcode that
    appeared in the analysis accession fragment file or RNA counts matrix. Not every barcode is
    guaranteed to be a real annotated cell, however. The QC columns are:
      - `analysis_set_accession` - the analysis accession the barcode came from
      - `barcode_sample` - the barcode
      - `annotated` - whether or not this barcode was annotated; may NOT be true for all barcode in this file
      - `found_in_rna` - whether or not this barcode was found in the RNA counts matrix
      - `found_in_atac` - whether or not this barcode was found in the ATAC fragments file
      - `pseudobulk_id` - either empty, or the pseudobulk that this barcode was assigned to
      - `rna_read_count` - the number of RNA reads
      - `gene_count` - the number of unique genes
      - `pct_mito` - the percent of reads from mitochrondrial genes
      - `pct_ribo` - the percent of reads from ribosomal genes
      - `num_frags` - the number of ATAC fragments
      - `pct_duplicated_reads` - the percent of ATAC reads that are duplicates
      - `nucleosomal_signal` - the fraction of mono-nucleosomal fragments (148-294 bp) to the fraction
        of nucleosome-free fragments (1-147 bp)
      - `tss_enrichment` - the enrichment of shifted insertions near TSS centers over the average
        number of shifted insertions at TSS flanks
  - `pseudobulk_qc.tsv.gz` - this file contains per-pseudobulk QC metrics. The QC metrics are not the
    mean or median of the corresponding QC metrics from the individual cells in the pseudobulk.
    Rather, these QC metrics are computed when treating the cluster as a single pseudobulked entity.
    The columns in this file are:
      - `directory_name` - the pseudobulk ID assigned by the pipeline
      - `pseudobulk` - a name obtained by replacing annotation_${annotation_index} with the `cell_name`.
      - `cell_name` - the value of `cell_name` as provided in the input annotations file.
      - `subsample` - the subsample.
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
- `analysis_accession_qc_reports/` — this directory will contain one QC file per analysis accession
  named `${analysis_accession}_per_cell_qc.tsv`. This QC file contains one row per barcode that
  appeared in the analysis accession fragment file or RNA counts matrix. Not every barcode is
  guaranteed to be a real annotated cell, however. The QC columns are:
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
    - `nucleosomal_signal` - the fraction of mono-nucleosomal fragments (148-294 bp) to the fraction
      of nucleosome-free fragments (1-147 bp)
    - `tss_enrichment` - the enrichment of shifted insertions near TSS centers over the average
      number of shifted insertions at TSS flanks
- `pseudobulk_qc.tsv` - this file contains per-pseudobulk QC metrics. The QC metrics are not the
  mean or median of the corresponding QC metrics from the individual cells in the pseudobulk.
  Rather, these QC metrics are computed when treating the cluster as a single pseudobulked entity.
  The columns in this file are:
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

## Annotation File Requirements

Annotations are required to adhere to the [single cell cell annotation table file format specification template](https://data.igvf.org/documents/a5f7c1d4-64a7-40fb-be0e-4b7d5596df01/)

## Genome Data

The `genome_data/` directory holds reference data for supported species (`human`, `mouse`). It
includes chromosome sizes, blacklist regions, gene information, and TSS data used during processing.

For details on how these files were generated, see `genome_data/README.md`.

## Processing and calculation notes

- Fragments that do not come from a chromosome in the species chromosome sizes file are
  **ignored entirely**. They will not even contribute towards the `num_frags` QC metric.
- The `tss_enrichment` is computed per barcode via the following procedure:
  - For each barcode, a generic window of 2001 bases (center +- 1000 bases) around TSS centers is
    instantiated. Index 0 corresponds to -1000 bases from all TSSs and index 2000 corresponds to
    +1000 bases from all TSSs.
  - From each ATAC fragment, two shifted (+4/-4) insertions are generated.
  - For each insertion, its distance to the closest TSS (based on the species TSS file) is computed.
    (chrM TSSs are ignored)
  - If the distance is within 1000 bp, the signed distance is computed (strand is taken into
    account). The value in the barcode's generic TSS window corresponding to that signed distance is
    incremented, representing finding an insertion at that distance from a TSS.
  - After all fragments have been processed, the average number of insertions per base in the
    flanking 200 bases of the generic TSS window are computed as
    `tss_insertions_flank_mean = (np.sum(tss_insertions[:100]) + np.sum(tss_insertions[-100:]))/200`.
  - Then, the average number of insertions per base in the central 5 bases are computed as
    `tss_insertions_center = np.mean(tss_insertions[1000-2:1000+3])`.
  - The TSS Enrichment is computed as the ratio between these values:
    `(tss_insertions_center/(tss_insertions_flank_mean if tss_insertions_flank_mean > 0 else 1))`.
- When computing `nucleosomal_signal`, the fragment size is computed based on the +4/-4 shifted
  insertion fragment endpoints rather than the raw fragment endpoints. Also, a pseudocount of 1 is
  applied to both mono-nucleosomal fragments (148-294 bp) and nucleosome-free fragments (1-147 bp).

## Developing

1. Changes should be made in your own `git` branch so that they can be merged via a pull request.
2. At the moment, minimal style guidelines are enforced with automated scripts and testing. Run
    ```bash
    pixi run checks
    ```
    and fix any problems that are found.

3. To test changes, you'll need to build docker images before running the pipeline. Run
    ```bash
    pixi run build-dockers
    ```
    Note that at the moment, this will clobber the dockers for everyone. A build system that's
    friendlier for multiple users will be coming soon.
