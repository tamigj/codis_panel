#!/bin/bash

#SBATCH --job-name=submit_jobs_0.75
#SBATCH --time=03:00:00
#SBATCH --mem=1G
#SBATCH --output=logs_0.75/submit_jobs_0.75.log
#SBATCH --error=logs_0.75/submit_jobs_0.75.err
#SBATCH -p hns,owners,normal
#SBATCH -c 1

source /scratch/groups/noahr/tami/codis_panel/config.sh

#---------------------------#
# EDIT FOR EACH EXPERIMENT  #
#---------------------------#
experiment='combinations'
n_snp_reps=10
n_ind_reps=10
fraction=0.75

ml R/4.2

mkdir -p logs_${fraction}
mkdir -p submit_${fraction}

# Define the specific snplist_ids and their parameters
declare -A snplist_params
snplist_params["0.05_NA_125000_NA"]="0.05 NA 125000 NA"
snplist_params["0.1_NA_125000_NA"]="0.1 NA 125000 NA"
snplist_params["NA_0_62500_NA"]="NA 0 62500 NA"
snplist_params["NA_0_125000_NA"]="NA 0 125000 NA"
snplist_params["NA_0.01_62500_NA"]="NA 0.01 62500 NA"
snplist_params["NA_0.05_125000_NA"]="NA 0.05 125000 NA"

# Loop through each snplist_id
for snplist_id in "${!snplist_params[@]}"; do

    # Parse parameters for this snplist_id
    read -r maf popmaf distance d <<< "${snplist_params[$snplist_id]}"

    for n_snps in $ALL_PANEL_SIZES_UP_TO_100;
    do

      for n_snp_rep in $(seq 1 $n_snp_reps);
      do

        # Make string for combo and snplist_id
        combo=${maf}_${popmaf}_${distance}_${d}
        snplist_id=${experiment}_${combo}_${n_snps}_${n_snp_rep}

        # 1. Prepare SNP list
        Rscript ${DIR_ANALYSES}/combinations/make_snplist_combos.R \
                ${experiment} ${n_snps} ${n_snp_rep} \
                ${maf} ${popmaf} ${distance} ${d}

        # 2. Submit the record matching pipeline
        for n_ind_rep in $(seq 1 $n_ind_reps);
        do

          n_run=${snplist_id}_${n_ind_rep}

          # Submit the record_matching_pipeline
          sbatch_script="submit_${n_run}.sh"

            # If results don't already exist
            if [[ ! -f "$DIR_OUTPUT_EXPERIMENTS/$fraction/run_$n_run/match_accuracies.csv" ]]; then

              # Create the sbatch script
              cat <<EOT > $sbatch_script
#!/bin/bash
#SBATCH --job-name=rm_${n_run}
#SBATCH --time=06:00:00
#SBATCH --mem=2G
#SBATCH --output=logs_${fraction}/rm_${n_run}.log
#SBATCH --error=logs_${fraction}/rm_${n_run}.err
#SBATCH -p hns,owners,normal
#SBATCH -c 1

$DIR_CODE/record_matching_pipeline_combinations.sh \
$snplist_id $n_ind_rep $fraction

EOT

              sbatch $sbatch_script
            fi

        done
      done
    done

done && mv submit*sh submit_$fraction/
