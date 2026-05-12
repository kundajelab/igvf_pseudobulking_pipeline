#!/usr/bin/env python3
"""
Generate IGVF portal submission TSVs for pseudobulk data.

Submission order:
  1. Document (QC files)
  2. Pseudobulk Set (links to QC docs)
  3. Tabular File, Matrix File, Signal File (link to pseudobulk sets)

Step 1: Generate pseudobulk_set and document TSVs
  python generate_igvf_submission.py --step 1 \
    --basedir /oak/stanford/groups/kasowski/sbaek1/igvf/igvf0/pseudobulks \
    --annotations /path/to/annotations_lookup.csv \
    --dataset IGVF0 \
    --input-file-sets IGVFDS4054XZLF,IGVFDS8289FYPK,IGVFDS1227ECUY \
    --output-name igvf0

Step 2: Generate file TSVs (run after step 1 objects are submitted)
  python generate_igvf_submission.py --step 2 \
    --basedir /oak/stanford/groups/kasowski/sbaek1/igvf/igvf0/pseudobulks \
    --dataset IGVF0 \
    --output-name igvf0 \
    --threads 8

  Files reference pseudobulk sets by alias, so no accession map is needed.
  Bed and tsv files marked needs_gzip are automatically gzipped (using pigz
  if available) before md5sum computation. The gzipped file is placed
  alongside the original. If only the .gz already exists, it is used as-is.
"""

import os
import hashlib
import argparse
import csv
import subprocess
import shutil

# ── Config ──────────────────────────────────────────────────────────────────
LAB = "/labs/anshul-kundaje/"
AWARD = "/awards/HG012069/"
FILE_SET_TYPE = "pseudobulk analysis"
ALIAS_PREFIX = "anshul-kundaje"
ASSEMBLY = "GRCh38"

# Reference files (genome + gene annotation)
REFERENCE_FILES = "IGVFFI0653VCGH,IGVFFI9573KOZR"

# ── File definitions ────────────────────────────────────────────────────────
TABULAR_FILES = {
    "fragments.tsv.gz": {
        "file_format": "tsv",
        "content_type": "fragments",
    },
    "peaks.narrowPeak": {
        "file_format": "bed",
        "file_format_type": "bed6+",
        "content_type": "peaks",
        "needs_gzip": True,
    },
    "pseudobulk_expression.tsv": {
        "file_format": "tsv",
        "content_type": "gene quantifications",
        "needs_gzip": True,
    },
}

MATRIX_FILES = {
    "rna_counts_mtx.h5ad": {
        "file_format": "h5ad",
        "content_type": "sparse gene count matrix",
        "principal_dimension": "cell",
        "secondary_dimensions": "gene",
    },
}

SIGNAL_FILES = {
    "peaks_minuslog10pval.bw": {
        "file_format": "bigWig",
        "content_type": "signal p-value",
        "strand_specificity": "unstranded",
    },
    "raw_insertions.bw": {
        "file_format": "bigWig",
        "content_type": "signal",
        "strand_specificity": "unstranded",
    },
}

DOCUMENT_FILES = {
    "per_cell_qc.tsv": {
        "document_type": "quality control report",
        "description": "Per-cell QC metrics for pseudobulk",
    },
}


def md5sum(filepath):
    h = hashlib.md5()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def gzip_file(filepath, threads=1):
    """Gzip a file in place using pigz (multi-core) or gzip. Returns path to .gz file."""
    gz_path = filepath + ".gz"
    if os.path.exists(gz_path):
        print(f"    Already gzipped: {gz_path}")
        return gz_path

    if shutil.which("pigz"):
        cmd = ["pigz", "-k", "-p", str(threads), filepath]
        print(f"    Gzipping with pigz ({threads} threads): {filepath}")
    else:
        cmd = ["gzip", "-k", filepath]
        print(f"    Gzipping with gzip: {filepath}")

    subprocess.run(cmd, check=True)
    if not os.path.exists(gz_path):
        raise RuntimeError(f"Gzip failed: {gz_path} not created")
    return gz_path


def resolve_gzip_upload(filepath, threads=1):
    """Return the .gz path to upload for a needs_gzip file.

    Handles three cases:
      - .gz already exists  -> use it (no re-gzip)
      - only original exists -> gzip it (keeps original)
      - neither exists       -> return None (caller should warn + skip)
    """
    gz_path = filepath + ".gz"
    if os.path.exists(gz_path):
        if not os.path.exists(filepath):
            print(f"    Already gzipped (original removed): {gz_path}")
        else:
            print(f"    Already gzipped: {gz_path}")
        return gz_path
    if os.path.exists(filepath):
        return gzip_file(filepath, threads=threads)
    return None


def load_annotations(csv_path, dataset):
    lookup = {}
    # utf-8-sig strips a leading BOM (Excel/Sheets exports); newline="" lets the
    # csv module handle CRLF line endings correctly.
    with open(csv_path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        reader.fieldnames = [name.strip() for name in reader.fieldnames]
        if "dataset" not in reader.fieldnames:
            raise ValueError(
                f"No 'dataset' column in {csv_path}. Found columns: {reader.fieldnames}"
            )
        for row in reader:
            if row["dataset"].strip() == dataset:
                ct = row["lab_celltype"].strip()
                cl_id = row["CL_ID"].strip()  # Keep as-is (e.g. CL_0000222)
                qualifier = row.get("qualifier", "").strip()
                lookup[ct] = {"cl_id": cl_id, "qualifier": qualifier}
    return lookup


def parse_folder_name(folder):
    parts = folder.split("-")
    if parts[0] != "annotation":
        return None, None
    sample_id = parts[-1]
    cell_type = "-".join(parts[1:-1])
    return cell_type, sample_id


def make_alias(dataset_tag, cell_type, sample_id, suffix=""):
    base = f"{ALIAS_PREFIX}:{dataset_tag}-{cell_type}-{sample_id}"
    if suffix:
        base += f"-{suffix}"
    return base


def scan_pseudobulk_dirs(basedir):
    entries = []
    for folder in sorted(os.listdir(basedir)):
        cell_type, sample_id = parse_folder_name(folder)
        if cell_type and sample_id:
            entries.append((folder, cell_type, sample_id))
    return entries


def write_tsv(outfile, rows):
    """Write rows as TSV without quoting (important for JSON fields like attachment)."""
    if not rows:
        return
    all_keys = list(dict.fromkeys(k for row in rows for k in row.keys()))
    with open(outfile, "w") as f:
        f.write("\t".join(all_keys) + "\n")
        for row in rows:
            values = [str(row.get(k, "")) for k in all_keys]
            f.write("\t".join(values) + "\n")
    print(f"\nWrote {len(rows)} rows to {outfile}")


def report_celltype_match(entries, ct_lookup):
    """Report bidirectional match between folder cell types and lookup entries.

    Returns the set of folder cell types that have no lookup entry (these become
    MISSING_ rows downstream).
    """
    folder_cts = {ct for _, ct, _ in entries}
    lookup_cts = set(ct_lookup)

    only_in_folders = folder_cts - lookup_cts   # folders with no lookup entry
    only_in_lookup = lookup_cts - folder_cts    # lookup rows never used by a folder

    print("\n── Cell type match report ──")
    print(f"  Folders: {len(folder_cts)} distinct cell types | "
          f"Lookup: {len(lookup_cts)} entries")

    if only_in_folders:
        print(f"  ⚠️  In folders but NOT in lookup ({len(only_in_folders)}) "
              f"-> will become MISSING_ rows: {sorted(only_in_folders)}")
    if only_in_lookup:
        print(f"  ⚠️  In lookup but NOT in any folder ({len(only_in_lookup)}) "
              f"-> unused lookup entries: {sorted(only_in_lookup)}")
    if not only_in_folders and not only_in_lookup:
        print(f"  ✓ Perfect match: all {len(folder_cts)} cell types matched "
              f"in both directions.")

    return only_in_folders


def step1_pseudobulk_and_docs(basedir, annotations_csv, dataset, input_file_sets,
                              outdir=".", output_name=None):
    """Generate document TSV and pseudobulk_set TSV (with document links)."""
    ct_lookup = load_annotations(annotations_csv, dataset)
    print(f"Loaded {len(ct_lookup)} cell type annotations for {dataset}:")
    for ct, info in ct_lookup.items():
        print(f"  {ct} -> {info['cl_id']}  qualifier={info['qualifier']}")

    entries = scan_pseudobulk_dirs(basedir)
    dataset_tag = dataset.lower().replace("#", "")
    prefix = f"{output_name}_" if output_name else ""

    # Bidirectional reconciliation between folders and lookup
    report_celltype_match(entries, ct_lookup)

    missing_ct = set()
    document_rows = []
    pseudobulk_rows = []

    for folder, cell_type, sample_id in entries:
        folder_path = os.path.join(basedir, folder)
        ps_alias = make_alias(dataset_tag, cell_type, sample_id)

        # --- Document (QC) ---
        doc_alias = ""
        for filename, fdef in DOCUMENT_FILES.items():
            filepath = os.path.join(folder_path, filename)
            if not os.path.exists(filepath):
                print(f"  WARNING: missing {filepath}")
                continue
            doc_alias = make_alias(dataset_tag, cell_type, sample_id,
                                   suffix=filename.replace(".", "_"))
            attachment = '{"path": "' + filepath + '"}'
            document_rows.append({
                "aliases": doc_alias,
                "award": AWARD,
                "lab": LAB,
                "document_type": fdef["document_type"],
                "description": f"{fdef['description']} for {cell_type} in {sample_id}",
                "attachment": attachment,
            })

        # --- Pseudobulk Set ---
        if cell_type in ct_lookup:
            cl_id = f"/sample-terms/{ct_lookup[cell_type]['cl_id']}/"
            qualifier = ct_lookup[cell_type]["qualifier"]
        else:
            cl_id = f"MISSING_{cell_type}"
            qualifier = ""
            missing_ct.add(cell_type)

        pseudobulk_rows.append({
            "aliases": ps_alias,
            "award": AWARD,
            "lab": LAB,
            "file_set_type": FILE_SET_TYPE,
            "cell_type": cl_id,
            "cell_qualifier": qualifier,
            "samples": sample_id,
            "input_file_sets": input_file_sets,
            "documents": doc_alias,
        })

    # Write document TSV
    if document_rows:
        out = os.path.join(outdir, f"{prefix}document.txt")
        write_tsv(out, document_rows)
        print(f"  Submit FIRST:")
        print(f"  iu_register -m staging -p document -i {out}")

    # Write pseudobulk_set TSV
    if pseudobulk_rows:
        out = os.path.join(outdir, f"{prefix}pseudobulk_set.txt")
        write_tsv(out, pseudobulk_rows)
        if missing_ct:
            print(f"\n⚠️  Missing CL IDs for cell types: {missing_ct}")
        print(f"  Submit SECOND:")
        print(f"  iu_register -m staging -p pseudobulk_set -i {out}")


def step2_file_tsvs(basedir, dataset, compute_md5=True,
                    outdir=".", output_name=None, threads=1):
    entries = scan_pseudobulk_dirs(basedir)
    dataset_tag = dataset.lower().replace("#", "")
    prefix = f"{output_name}_" if output_name else ""

    tabular_rows = []
    matrix_rows = []
    signal_rows = []

    for folder, cell_type, sample_id in entries:
        file_set = make_alias(dataset_tag, cell_type, sample_id)
        folder_path = os.path.join(basedir, folder)

        # --- Tabular files (assembly + controlled_access required) ---
        for filename, fdef in TABULAR_FILES.items():
            filepath = os.path.join(folder_path, filename)

            if fdef.get("needs_gzip"):
                upload_path = resolve_gzip_upload(filepath, threads=threads)
                if upload_path is None:
                    print(f"  WARNING: missing {filepath} (and {filepath}.gz)")
                    continue
            else:
                if not os.path.exists(filepath):
                    print(f"  WARNING: missing {filepath}")
                    continue
                upload_path = filepath

            file_alias = make_alias(dataset_tag, cell_type, sample_id,
                                    suffix=filename.replace(".", "_"))
            md5 = ""
            if compute_md5:
                print(f"  Computing md5: {upload_path}")
                md5 = md5sum(upload_path)
            row = {
                "aliases": file_alias,
                "award": AWARD,
                "lab": LAB,
                "file_set": file_set,
                "file_format": fdef["file_format"],
                "content_type": fdef["content_type"],
                "controlled_access": "false",
                # "assembly": ASSEMBLY,
                "md5sum": md5,
                "submitted_file_name": upload_path,
                "reference_files": REFERENCE_FILES,
            }
            if "file_format_type" in fdef:
                row["file_format_type"] = fdef["file_format_type"]
            tabular_rows.append(row)

        # --- Matrix files (no assembly, no controlled_access) ---
        for filename, fdef in MATRIX_FILES.items():
            filepath = os.path.join(folder_path, filename)
            if not os.path.exists(filepath):
                print(f"  WARNING: missing {filepath}")
                continue
            file_alias = make_alias(dataset_tag, cell_type, sample_id,
                                    suffix=filename.replace(".", "_"))
            md5 = ""
            if compute_md5:
                print(f"  Computing md5: {filepath}")
                md5 = md5sum(filepath)
            matrix_rows.append({
                "aliases": file_alias,
                "award": AWARD,
                "lab": LAB,
                "file_set": file_set,
                "file_format": fdef["file_format"],
                "content_type": fdef["content_type"],
                "md5sum": md5,
                "submitted_file_name": filepath,
                "reference_files": REFERENCE_FILES,
                "principal_dimension": fdef["principal_dimension"],
                "secondary_dimensions": fdef["secondary_dimensions"],
            })

        # --- Signal files (no assembly, no controlled_access, strand_specificity required) ---
        for filename, fdef in SIGNAL_FILES.items():
            filepath = os.path.join(folder_path, filename)
            if not os.path.exists(filepath):
                print(f"  WARNING: missing {filepath}")
                continue
            file_alias = make_alias(dataset_tag, cell_type, sample_id,
                                    suffix=filename.replace(".", "_"))
            md5 = ""
            if compute_md5:
                print(f"  Computing md5: {filepath}")
                md5 = md5sum(filepath)
            signal_rows.append({
                "aliases": file_alias,
                "award": AWARD,
                "lab": LAB,
                "file_set": file_set,
                "file_format": fdef["file_format"],
                "content_type": fdef["content_type"],
                "strand_specificity": fdef["strand_specificity"],
                "md5sum": md5,
                "submitted_file_name": filepath,
                "reference_files": REFERENCE_FILES,
            })

    if tabular_rows:
        out = os.path.join(outdir, f"{prefix}tabular_file.txt")
        write_tsv(out, tabular_rows)
        print(f"  iu_register -m staging -p tabular_file -i {out}")

    if matrix_rows:
        out = os.path.join(outdir, f"{prefix}matrix_file.txt")
        write_tsv(out, matrix_rows)
        print(f"  iu_register -m staging -p matrix_file -i {out}")

    if signal_rows:
        out = os.path.join(outdir, f"{prefix}signal_file.txt")
        write_tsv(out, signal_rows)
        print(f"  iu_register -m staging -p signal_file -i {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--step", type=int, required=True, choices=[1, 2])
    parser.add_argument("--basedir", required=True)
    parser.add_argument("--annotations", help="Path to annotations CSV (step 1)")
    parser.add_argument("--dataset", required=True, help="Dataset name, e.g. IGVF0")
    parser.add_argument("--input-file-sets", help="Comma-separated analysis set IDs (step 1)")
    parser.add_argument("--no-md5", action="store_true", help="Skip md5 computation")
    parser.add_argument("--outdir", default=".", help="Output directory for TSVs")
    parser.add_argument("--output-name", help="Step 1: prefix for document + pseudobulk_set files. "
                        "Step 2: prefix for file TSVs (e.g. igvf0 -> igvf0_tabular_file.txt)")
    parser.add_argument("--threads", type=int, default=1,
                        help="Number of threads for pigz compression (default: 1)")
    args = parser.parse_args()

    if args.step == 1:
        if not args.annotations:
            parser.error("--annotations required for step 1")
        ifs = args.input_file_sets or ""
        step1_pseudobulk_and_docs(args.basedir, args.annotations, args.dataset, ifs,
                                  args.outdir, args.output_name)
    elif args.step == 2:
        step2_file_tsvs(args.basedir, args.dataset,
                        compute_md5=not args.no_md5, outdir=args.outdir,
                        output_name=args.output_name, threads=args.threads)