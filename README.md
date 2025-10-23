# IGVF Pseudobulking (Preliminary Repository)

This repository is a **PRELIMINARY** tool to assist with IGVF pseudobulking.  
It helps you download IGVF single-cell data and perform pseudobulking.

---

## Instructions

### 0) Prerequisites
- Ensure you have **download access** to data from the IGVF portal.  
- You’ll need your **access key** and **secret key**.

### 1) Clone this repository
```bash
git clone https://github.com/yourusername/igvf-pseudobulk.git
cd igvf-pseudobulk
```

### 2) Install the environment
Use the provided [environment.yml](environment.yml) file to create the conda environment.
```bash
conda env create -f environment.yml
conda activate igvf-pseudobulk
```

### 3) Choose your workspace directory
```bash
workspace=/path/to/dir
```

### 4) Download IGVF data
To download IGVF data, first identify **analysis accessions** that you want to process together.  
These should correspond to analyses that you will pseudobulk together.

For example, suppose we want to work with the **IGVF1** dataset (lab accession `IGVF1`).  
This corresponds to the uniform analysis accessions `IGVFDS6430MYNQ` and `IGVFDS5417HJRJ`.

Run the following:
```bash
bash igvf_process.sh IGVFDS6430MYNQ,IGVFDS5417HJRJ ${access_key} ${secret_key} ${workspace}
```

This will download the RNA count matrices and fragment files corresponding to the scRNA and scATAC data.

> From here on, you can do your own pseudobulking. We **strongly encourage** using this standard pipeline code, but any pseudobulking approach can be used.

---

### 5) Prepare metadata for pseudobulking
Next, prepare a **single metadata file** for both analyses.

The metadata file must include:
- a column called `barcode` → ACGT string corresponding to a cell
- a column called `analysis_accession` → one of the accessions

Each `(barcode, analysis_accession)` pair specifies a single cell.

You can include any number of **additional columns** that define pseudobulk groupings (e.g., by cell type, time point, etc.).

See [`test_metadata.tsv`](test_metadata.tsv) for an example metadata file.

> **NOTE:** If you have one annotation file per analysis accession, consider doing the following:
> 1. Ensure each file has a `barcode` column.  
> 2. Add a column `uniform_accession` with the same value for all rows in that file.  
> 3. Merge them into a single annotation file.

---

### 6) Pseudobulk processing (preliminary)

> **NOTE: This part is still preliminary and not yet an IGVF official standard.**  
> **RNA pseudobulking is not yet implemented (as of 10/23/2025).**

In this step, fragments will be pseudobulked and peaks will be called.

You need the following files (paths are from the repo root):
```bash
chr_order_file="chr_info_data/GRCh38_EBV_sorted_standard.chrom.sizes.tsv"
blacklist_file="chr_info_data/hg38.blacklist.bed.gz"
```

Then, determine how many CPUs to use for parallelization, for example:
```bash
numcpus=4
```

Run the pseudobulk pipeline:
```bash
bash igvf_process.sh ${workspace} ${metadata_path} ${chr_order_file} ${blacklist_file} ${numcpus}
```

This process will take several minutes to tens of minutes.  
Upon completion, you will have the following directory structure:

```bash
{workspace}/pseudobulked_fragments/   # contains fragments as .tsv files
{workspace}/peaks/                    # contains narrowPeak peaks + MACS3 signal tracks as bigWigs
```

---

