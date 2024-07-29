#!/bin/bash

source /scratch/groups/noahr/tami/codis_panel/config.sh

#---------------------------#
# EDIT FOR EACH EXPERIMENT  #
#---------------------------#
experiment='random'
n_snps_per_str=$ALL_PANEL_SIZES_CSV
n_snp_reps=10
n_ind_reps=10
fraction=0.25

###############################################################
ml R/4.2

Rscript $DIR_CODE/summarize_results.R \
    $CODIS_STRS_RM \
    $experiment \
    $n_snps_per_str $n_snp_reps $n_ind_reps $fraction \
    $DIR_DATA_SNP_LISTS $DIR_OUTPUT_EXPERIMENTS $DIR_OUTPUT_RM_SUMMARIES
