

# IGVF0 
python /oak/stanford/groups/akundaje/sbaek1/command/IGVF2025/portal_upload/generate_igvf_submission_file_20260407.py --step 1 \
  --basedir /oak/stanford/groups/kasowski/sbaek1/igvf/igvf0/pseudobulks \
  --annotations "/oak/stanford/groups/kasowski/sbaek1/igvf/lab_annotation_lookup.csv" \
  --dataset IGVF0 \
  --input-file-sets IGVFDS1612ZNCA \
  --outdir /oak/stanford/groups/kasowski/sbaek1/igvf/submission \
  --output-name igvf0

python /oak/stanford/groups/akundaje/sbaek1/command/IGVF2025/portal_upload/generate_igvf_submission_file_20260407.py --step 2 \
  --basedir /oak/stanford/groups/kasowski/sbaek1/igvf/igvf0/pseudobulks \
  --dataset IGVF0 \
  --outdir /oak/stanford/groups/kasowski/sbaek1/igvf/submission \
  --output-name igvf0 \
  --threads 20


python /oak/stanford/groups/akundaje/sbaek1/command/IGVF2025/portal_upload/generate_igvf_submission_file_20260407.py --step 1 \
  --basedir /oak/stanford/groups/kasowski/sbaek1/igvf/igvf1/pseudobulks \
  --annotations "/oak/stanford/groups/kasowski/sbaek1/igvf/lab_annotation_lookup_20260417.csv" \
  --dataset IGVF1 \
  --input-file-sets IGVFDS7786JSWZ \
  --outdir /oak/stanford/groups/kasowski/sbaek1/igvf/submission \
  --output-name igvf1

python /oak/stanford/groups/akundaje/sbaek1/command/IGVF2025/portal_upload/generate_igvf_submission_file_20260407.py --step 2 \
  --basedir /oak/stanford/groups/kasowski/sbaek1/igvf/igvf1/pseudobulks \
  --dataset IGVF1 \
  --outdir /oak/stanford/groups/kasowski/sbaek1/igvf/submission \
  --output-name igvf1 \
  --threads 20

### IGVF 1
python /oak/stanford/groups/akundaje/sbaek1/command/IGVF2025/portal_upload/generate_igvf_submission_file_20260407.py --step 1 \
  --basedir /oak/stanford/groups/kasowski/sbaek1/igvf/igvf1/pseudobulks \
  --annotations "/oak/stanford/groups/kasowski/sbaek1/igvf/lab_annotation_lookup_20260507.csv" \
  --dataset IGVF1 \
  --input-file-sets IGVFDS7786JSWZ \
  --outdir /oak/stanford/groups/kasowski/sbaek1/igvf/submission \
  --output-name igvf1

python /oak/stanford/groups/akundaje/sbaek1/command/IGVF2025/portal_upload/generate_igvf_submission_file_20260407.py --step 2 \
  --basedir /oak/stanford/groups/kasowski/sbaek1/igvf/igvf1/pseudobulks \
  --dataset IGVF1 \
  --outdir /oak/stanford/groups/kasowski/sbaek1/igvf/submission \
  --output-name igvf1 \
  --threads 20

### IGVF 3
python /oak/stanford/groups/akundaje/sbaek1/command/IGVF2025/portal_upload/generate_igvf_submission_file_20260407.py --step 1 \
  --basedir /oak/stanford/groups/kasowski/sbaek1/igvf/igvf3/pseudobulks \
  --annotations "/oak/stanford/groups/kasowski/sbaek1/igvf/lab_annotation_lookup_20260507.csv" \
  --dataset IGVF3 \
  --input-file-sets IGVFDS1170FKEK \
  --outdir /oak/stanford/groups/kasowski/sbaek1/igvf/submission \
  --output-name igvf3

python /oak/stanford/groups/akundaje/sbaek1/command/IGVF2025/portal_upload/generate_igvf_submission_file_20260407.py --step 2 \
  --basedir /oak/stanford/groups/kasowski/sbaek1/igvf/igvf3/pseudobulks \
  --dataset IGVF3 \
  --outdir /oak/stanford/groups/kasowski/sbaek1/igvf/submission \
  --output-name igvf3 \
  --threads 20



### 20260603 IGVF2,4,5,6,10,12,17
# IGVF2
python /oak/stanford/groups/akundaje/sbaek1/command/IGVF2025/portal_upload/generate_igvf_submission_file_20260407.py --step 1 \
  --basedir /oak/stanford/groups/kasowski/sbaek1/igvf/igvf2/pseudobulks \
  --annotations "/oak/stanford/groups/kasowski/sbaek1/igvf/lab_annotation_lookup_20260603.csv" \
  --dataset IGVF2 \
  --input-file-sets IGVFDS4258QDIT \
  --outdir /oak/stanford/groups/kasowski/sbaek1/igvf/submission \
  --output-name igvf2

python /oak/stanford/groups/akundaje/sbaek1/command/IGVF2025/portal_upload/generate_igvf_submission_file_20260407.py --step 2 \
  --basedir /oak/stanford/groups/kasowski/sbaek1/igvf/igvf2/pseudobulks \
  --dataset IGVF2 \
  --outdir /oak/stanford/groups/kasowski/sbaek1/igvf/submission \
  --output-name igvf2 \
  --threads 10

# IGVF4
python /oak/stanford/groups/akundaje/sbaek1/command/IGVF2025/portal_upload/generate_igvf_submission_file_20260407.py --step 1 \
  --basedir /oak/stanford/groups/kasowski/sbaek1/igvf/igvf4/pseudobulks \
  --annotations "/oak/stanford/groups/kasowski/sbaek1/igvf/lab_annotation_lookup_20260603.csv" \
  --dataset IGVF4 \
  --input-file-sets IGVFDS5875AFXS \
  --outdir /oak/stanford/groups/kasowski/sbaek1/igvf/submission \
  --output-name igvf4

python /oak/stanford/groups/akundaje/sbaek1/command/IGVF2025/portal_upload/generate_igvf_submission_file_20260407.py --step 2 \
  --basedir /oak/stanford/groups/kasowski/sbaek1/igvf/igvf4/pseudobulks \
  --dataset IGVF4 \
  --outdir /oak/stanford/groups/kasowski/sbaek1/igvf/submission \
  --output-name igvf4 \
  --threads 10

# IGVF5
python /oak/stanford/groups/akundaje/sbaek1/command/IGVF2025/portal_upload/generate_igvf_submission_file_20260407.py --step 1 \
  --basedir /oak/stanford/groups/kasowski/sbaek1/igvf/igvf5/pseudobulks \
  --annotations "/oak/stanford/groups/kasowski/sbaek1/igvf/lab_annotation_lookup_20260603.csv" \
  --dataset IGVF5 \
  --input-file-sets IGVFDS5851LHUQ \
  --outdir /oak/stanford/groups/kasowski/sbaek1/igvf/submission \
  --output-name igvf5

python /oak/stanford/groups/akundaje/sbaek1/command/IGVF2025/portal_upload/generate_igvf_submission_file_20260407.py --step 2 \
  --basedir /oak/stanford/groups/kasowski/sbaek1/igvf/igvf5/pseudobulks \
  --dataset IGVF5 \
  --outdir /oak/stanford/groups/kasowski/sbaek1/igvf/submission \
  --output-name igvf5 \
  --threads 10

# IGVF6
python /oak/stanford/groups/akundaje/sbaek1/command/IGVF2025/portal_upload/generate_igvf_submission_file_20260407.py --step 1 \
  --basedir /oak/stanford/groups/kasowski/sbaek1/igvf/igvf6/pseudobulks \
  --annotations "/oak/stanford/groups/kasowski/sbaek1/igvf/lab_annotation_lookup_20260603.csv" \
  --dataset IGVF6 \
  --input-file-sets IGVFDS4000IHRJ \
  --outdir /oak/stanford/groups/kasowski/sbaek1/igvf/submission \
  --output-name igvf6

python /oak/stanford/groups/akundaje/sbaek1/command/IGVF2025/portal_upload/generate_igvf_submission_file_20260407.py --step 2 \
  --basedir /oak/stanford/groups/kasowski/sbaek1/igvf/igvf6/pseudobulks \
  --dataset IGVF6 \
  --outdir /oak/stanford/groups/kasowski/sbaek1/igvf/submission \
  --output-name igvf6 \
  --threads 10

# IGVF10
python /oak/stanford/groups/akundaje/sbaek1/command/IGVF2025/portal_upload/generate_igvf_submission_file_20260407.py --step 1 \
  --basedir /oak/stanford/groups/kasowski/sbaek1/igvf/igvf10/pseudobulks \
  --annotations "/oak/stanford/groups/kasowski/sbaek1/igvf/lab_annotation_lookup_20260603.csv" \
  --dataset IGVF10 \
  --input-file-sets IGVFDS3549FUQU \
  --outdir /oak/stanford/groups/kasowski/sbaek1/igvf/submission \
  --output-name igvf10

python /oak/stanford/groups/akundaje/sbaek1/command/IGVF2025/portal_upload/generate_igvf_submission_file_20260407.py --step 2 \
  --basedir /oak/stanford/groups/kasowski/sbaek1/igvf/igvf10/pseudobulks \
  --dataset IGVF10 \
  --outdir /oak/stanford/groups/kasowski/sbaek1/igvf/submission \
  --output-name igvf10 \
  --threads 10


# IGVF12
python /oak/stanford/groups/akundaje/sbaek1/command/IGVF2025/portal_upload/generate_igvf_submission_file_20260407.py --step 1 \
  --basedir /oak/stanford/groups/kasowski/sbaek1/igvf/igvf12/pseudobulks \
  --annotations "/oak/stanford/groups/kasowski/sbaek1/igvf/lab_annotation_lookup_20260603.csv" \
  --dataset IGVF12 \
  --input-file-sets IGVFDS8692JVTD \
  --outdir /oak/stanford/groups/kasowski/sbaek1/igvf/submission \
  --output-name igvf12

python /oak/stanford/groups/akundaje/sbaek1/command/IGVF2025/portal_upload/generate_igvf_submission_file_20260407.py --step 2 \
  --basedir /oak/stanford/groups/kasowski/sbaek1/igvf/igvf12/pseudobulks \
  --dataset IGVF12 \
  --outdir /oak/stanford/groups/kasowski/sbaek1/igvf/submission \
  --output-name igvf12 \
  --threads 10

# IGVF17
python /oak/stanford/groups/akundaje/sbaek1/command/IGVF2025/portal_upload/generate_igvf_submission_file_20260407.py --step 1 \
  --basedir /oak/stanford/groups/kasowski/sbaek1/igvf/igvf17/pseudobulks \
  --annotations "/oak/stanford/groups/kasowski/sbaek1/igvf/lab_annotation_lookup_20260603.csv" \
  --dataset IGVF17 \
  --input-file-sets IGVFDS7039GRZR \
  --outdir /oak/stanford/groups/kasowski/sbaek1/igvf/submission \
  --output-name igvf17

python /oak/stanford/groups/akundaje/sbaek1/command/IGVF2025/portal_upload/generate_igvf_submission_file_20260407.py --step 2 \
  --basedir /oak/stanford/groups/kasowski/sbaek1/igvf/igvf17/pseudobulks \
  --dataset IGVF17 \
  --outdir /oak/stanford/groups/kasowski/sbaek1/igvf/submission \
  --output-name igvf17 \
  --threads 10  