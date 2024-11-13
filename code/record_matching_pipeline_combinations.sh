#!/bin/bash

snplist_id=$1
n_ind_rep=$2
fraction=$3

SECONDS=0

# Construct variables
n_run=${snplist_id}_${n_ind_rep}

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
for str in $CODIS_STRS;
do
  $DIR_CODE/B_impute_STR_from_SNPs.sh $DIR_DATA_NRUN \
                                      $DIR_OUTPUT_NRUN \
                                      $str
done &&

  # 3. Run record matching
ml R/4.2
Rscript $DIR_CODE/C_compute_match_accuracies.R $CODIS_STRS_RM \
                                               $BEAGLE_JAR \
                                               $VCF_EXE \
                                               $DIR_DATA_NRUN \
                                               $DIR_OUTPUT_NRUN

# Calculate hours, minutes, and seconds
hours=$((SECONDS / 3600))
minutes=$(( (SECONDS % 3600) / 60 ))
seconds=$((SECONDS % 60))

echo "Elapsed time: ${hours}h ${minutes}m ${seconds}s"
