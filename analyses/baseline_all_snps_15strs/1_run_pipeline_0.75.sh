#!/bin/bash

#SBATCH --job-name=submit_jobs_0.75
#SBATCH --time=00:30:00
#SBATCH --mem=1G
#SBATCH --output=logs_0.75/submit_jobs_0.75.log
#SBATCH --error=logs_0.75/submit_jobs_0.75.err
#SBATCH -p hns,owners,normal
#SBATCH -c 1

#---------------------------#
# EDIT FOR EACH EXPERIMENT  #
#---------------------------#
experiment='random_15str'
n_snps='all'
n_snp_reps=1
n_ind_reps=100
fraction=0.75

###############################################################
source /scratch/groups/noahr/tami/codis_panel/config.sh

ml R/4.2

mkdir -p logs_${fraction}
mkdir -p submit_${fraction}

for n_snp_rep in $(seq 1 $n_snp_reps);
do

  # 1. Prepare SNP list
  snplist_id=${experiment}_${n_snps}_${n_snp_rep}

  Rscript make_snplists_all_snps.R ${snplist_id}

  # 2. Submit the record matching pipeline
  for n_ind_rep in $(seq 1 $n_ind_reps);
  do

    n_run=${snplist_id}_${n_ind_rep}

    # Submit the record_matching_pipeline
    sbatch_script="submit_${n_run}.sh"

    # Create the sbatch script
    cat <<EOT > $sbatch_script
#!/bin/bash
#SBATCH --job-name=rm_${n_run}
#SBATCH --time=05:00:00
#SBATCH --mem=2G
#SBATCH --output=logs_${fraction}/rm_${n_run}.log
#SBATCH --error=logs_${fraction}/rm_${n_run}.err
#SBATCH -p hns,owners,normal
#SBATCH -c 1

$DIR_CODE/record_matching_pipeline_15strs.sh \
$experiment $n_snps $n_snp_rep $n_ind_rep $fraction

EOT

    # Submit the sbatch script
    sbatch $sbatch_script

  done
done && mv submit*sh submit_${fraction}/
