#!/bin/bash

source /scratch/groups/noahr/tami/codis_panel/config.sh

#---------------------------#
# EDIT FOR EACH EXPERIMENT  #
#---------------------------#
experiment='random_15str'
n_snps_per_str='all'
n_snp_reps=1
n_ind_reps=100
fraction=0.75

###############################################################
ml R/4.2

Rscript $DIR_CODE/summarize_results.R \
    $CODIS_15_STRS_RM \
    $experiment \
    $n_snps_per_str $n_snp_reps $n_ind_reps $fraction \
    $DIR_DATA_SNP_LISTS $DIR_OUTPUT_EXPERIMENTS $DIR_OUTPUT_RM_SUMMARIES
