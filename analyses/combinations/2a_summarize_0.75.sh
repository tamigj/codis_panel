#!/bin/bash

source /scratch/groups/noahr/tami/codis_panel/config.sh

#---------------------------#
# EDIT FOR EACH EXPERIMENT  #
#---------------------------#
experiment='combinations'
combos_file="$DIR_OUTPUT_SUMSTATS/feasible_combinations_2_100.csv"
n_snps_per_str=$ALL_PANEL_SIZES_UP_TO_100_CSV
n_snp_reps=10
n_ind_reps=10
fraction=0.75

###############################################################
ml R/4.2

Rscript $DIR_CODE/summarize_results_combinations.R \
    $CODIS_STRS_RM \
    $experiment $combos_file \
    $n_snps_per_str $n_snp_reps $n_ind_reps $fraction \
    $DIR_DATA_SNP_LISTS $DIR_OUTPUT_EXPERIMENTS $DIR_OUTPUT_RM_SUMMARIES
