conda activate igvf_utils_new
export IGVF_API_KEY=
export IGVF_SECRET_KEY=
export IGVF_LAB='anshul-kundaje'
export IGVF_AWARD='HG012069'



## staging
export IGVF_API_KEY=
export IGVF_SECRET_KEY=

cd /oak/stanford/groups/kasowski/sbaek1/igvf/submission/



### IGVF0 (test upload to staging)
iu_register -m staging -p document -i igvf0_document.txt          # 1st
iu_register -m staging -p pseudobulk_set -i igvf0_pseudobulk_set.txt  # 2nd
iu_register -m staging -p tabular_file -i igvf0_tabular_file.txt  # 3rd
iu_register -m staging -p matrix_file -i igvf0_matrix_file.txt    # 4th
iu_register -m staging -p signal_file -i igvf0_signal_file.txt    # 5th
iu_register --patch -m staging -p pseudobulk_set -w -i igvf0_pseudobulk_set_fixed.txt  # 6th

### IGVF1 
iu_register -m staging -p document -i igvf1_document.txt          # 1st
iu_register -m staging -p pseudobulk_set -i igvf1_pseudobulk_set.txt  # 2nd
iu_register -m staging -p tabular_file -i igvf1_tabular_file.txt  # 3rd
iu_register -m staging -p matrix_file -i igvf1_matrix_file.txt    # 4th
iu_register -m staging -p signal_file -i igvf1_signal_file.txt    # 5th

#### actual upload (manually removed barcode with 0)
## IGVF1
iu_register -m prod -p document -i igvf1_document.txt          # 1st
iu_register -m prod -p pseudobulk_set -i igvf1_pseudobulk_set.txt  # 2nd
cut -f1-7,9-12 igvf1_tabular_file.txt > igvf1_tabular_file_noassembly.txt ## remove assembly column
iu_register -m prod -p tabular_file -i igvf1_tabular_file_noassembly.txt  # 3rd
iu_register -m prod -p matrix_file -i igvf1_matrix_file.txt    # 4th
iu_register -m prod -p signal_file -i igvf1_signal_file.txt    # 5th

## IGVF3
iu_register -m prod -p document -i igvf3_document.txt          # 1st
iu_register -m prod -p pseudobulk_set -i igvf3_pseudobulk_set.txt  # 2nd
iu_register -m prod -p tabular_file -i igvf3_tabular_file.txt  # 3rd
iu_register -m prod -p matrix_file -i igvf3_matrix_file.txt    # 4th
iu_register -m prod -p signal_file -i igvf3_signal_file.txt    # 5th

### 20260603
## IGVF2
iu_register -m prod -p document -i igvf2_document.txt          # 1st
iu_register -m prod -p pseudobulk_set -i igvf2_pseudobulk_set.txt  # 2nd
iu_register -m prod -p tabular_file -i igvf2_tabular_file.txt  # 3rd
iu_register -m prod -p matrix_file -i igvf2_matrix_file.txt    # 4th
iu_register -m prod -p signal_file -i igvf2_signal_file.txt    # 5th

## IGVF4
iu_register -m prod -p document -i igvf4_document.txt          # 1st
iu_register -m prod -p pseudobulk_set -i igvf4_pseudobulk_set.txt  # 2nd
iu_register -m prod -p tabular_file -i igvf4_tabular_file.txt  # 3rd
iu_register -m prod -p matrix_file -i igvf4_matrix_file.txt    # 4th
iu_register -m prod -p signal_file -i igvf4_signal_file.txt    # 5th

## IGVF5
iu_register -m prod -p document -i igvf5_document.txt          # 1st
iu_register -m prod -p pseudobulk_set -i igvf5_pseudobulk_set.txt  # 2nd
iu_register -m prod -p tabular_file -i igvf5_tabular_file.txt  # 3rd
iu_register -m prod -p matrix_file -i igvf5_matrix_file.txt    # 4th
iu_register -m prod -p signal_file -i igvf5_signal_file.txt    # 5th

## IGVF6
iu_register -m prod -p document -i igvf6_document.txt          # 1st
iu_register -m prod -p pseudobulk_set -i igvf6_pseudobulk_set.txt  # 2nd
iu_register -m prod -p tabular_file -i igvf6_tabular_file.txt  # 3rd
iu_register -m prod -p matrix_file -i igvf6_matrix_file.txt    # 4th
iu_register -m prod -p signal_file -i igvf6_signal_file.txt    # 5th

## IGVF10
iu_register -m prod -p document -i igvf10_document.txt          # 1st
iu_register -m prod -p pseudobulk_set -i igvf10_pseudobulk_set.txt  # 2nd
iu_register -m prod -p tabular_file -i igvf10_tabular_file.txt  # 3rd
iu_register -m prod -p matrix_file -i igvf10_matrix_file.txt    # 4th
iu_register -m prod -p signal_file -i igvf10_signal_file.txt    # 5th

## IGVF12
iu_register -m prod -p document -i igvf12_document.txt          # 1st
iu_register -m prod -p pseudobulk_set -i igvf12_pseudobulk_set.txt  # 2nd
iu_register -m prod -p tabular_file -i igvf12_tabular_file.txt  # 3rd
iu_register -m prod -p matrix_file -i igvf12_matrix_file.txt    # 4th
iu_register -m prod -p signal_file -i igvf12_signal_file.txt    # 5th

## IGVF17
iu_register -m prod -p document -i igvf17_document.txt          # 1st
iu_register -m prod -p pseudobulk_set -i igvf17_pseudobulk_set.txt  # 2nd
iu_register -m prod -p tabular_file -i igvf17_tabular_file.txt  # 3rd
iu_register -m prod -p matrix_file -i igvf17_matrix_file.txt    # 4th
iu_register -m prod -p signal_file -i igvf17_signal_file.txt    # 5th