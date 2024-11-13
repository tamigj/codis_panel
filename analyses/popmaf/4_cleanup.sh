#!/bin/bash

#SBATCH --job-name=cleanup_popmaf
#SBATCH --time=01:00:00
#SBATCH --mem=1G
#SBATCH --output=cleanup_popmaf.log
#SBATCH --error=cleanup_popmaf.err
#SBATCH -p hns,owners,normal
#SBATCH -c 1

source /scratch/groups/noahr/tami/codis_panel/config.sh

#---------------------------#
# EDIT FOR EACH EXPERIMENT  #
#---------------------------#
experiment='popmaf'
n_snps_per_str=$SELECTED_PANEL_SIZES
n_snp_reps=10
n_ind_reps=10

###############################################################

for popmaf_filter in $POPMAF_FILTERS;
do

  for n_snps in $n_snps_per_str;
  do

    for n_snp_rep in $(seq 1 $n_snp_reps);
    do

      for n_ind_rep in $(seq 1 $n_ind_reps);
      do

        for fraction in $FRACTIONS;
        do

          snplist_id=${experiment}_${popmaf_filter}_${n_snps}_${n_snp_rep}
          n_run=${snplist_id}_${n_ind_rep}

          $DIR_CODE/cleanup_intermediate_files.sh \
          $snplist_id $n_run $fraction

        done
      done
    done
  done
done
