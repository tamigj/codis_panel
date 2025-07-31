#!/bin/bash

SECONDS=0

experiment=$1
n_snps=$2
n_snp_rep=$3
n_ind_rep=$4
fraction=$5

# Construct variables
snplist_id=${experiment}_${n_snps}_${n_snp_rep}
n_run=${snplist_id}_${n_ind_rep}

source /scratch/groups/noahr/tami/codis_panel/config.sh

# Update directories
DIR_DATA_PARTITIONS_FRACTIONS="$DIR_DATA_PARTITIONS/$fraction"
DIR_DATA_NRUN="$DIR_DATA_EXPERIMENTS/$fraction/run_$n_run"
DIR_OUTPUT_NRUN="$DIR_OUTPUT_EXPERIMENTS/$fraction/run_$n_run"


# 1. Prepare files
$DIR_CODE/A_prepare_files.sh $n_ind_rep \
                             $snplist_id \
                             $DIR_DATA_NRUN \
                             $DIR_DATA_PARTITIONS_FRACTIONS &&

# 2. Impute STRs
for str in $CODIS_15_STRS;
do
  $DIR_CODE/B_impute_STR_from_SNPs.sh $DIR_DATA_NRUN \
                                      $DIR_OUTPUT_NRUN \
                                      $str
done &&

  # 3. Run record matching
ml R/4.2
Rscript $DIR_CODE/C_compute_match_accuracies.R $CODIS_15_STRS_RM \
                                               $BEAGLE_JAR \
                                               $VCF_EXE \
                                               $DIR_DATA_NRUN \
                                               $DIR_OUTPUT_NRUN

# Calculate hours, minutes, and seconds
hours=$((SECONDS / 3600))
minutes=$(( (SECONDS % 3600) / 60 ))
seconds=$((SECONDS % 60))

echo "Elapsed time: ${hours}h ${minutes}m ${seconds}s"
