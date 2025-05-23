#!/bin/bash

#SBATCH --job-name=rerun_jobs_0.75
#SBATCH --time=03:00:00
#SBATCH --mem=1G
#SBATCH --output=logs_0.75/rerun_jobs_0.75.log
#SBATCH --error=logs_0.75/rerun_jobs_0.75.err
#SBATCH -p hns,owners,normal
#SBATCH -c 1

source /scratch/groups/noahr/tami/codis_panel/config.sh

#---------------------------#
# EDIT FOR EACH EXPERIMENT  #
#---------------------------#
experiment='combinations'
n_snps_per_str=$ALL_PANEL_SIZES_UP_TO_100
n_snp_reps=10
n_ind_reps=10
fraction=0.75

time_for_job="15:00:00"

# Directory paths
submit_dir="submit_$fraction"
out_dir="$DIR_OUTPUT_EXPERIMENTS/$fraction"

# Main loop
combos_file="$DIR_OUTPUT_SUMSTATS/feasible_combinations_2_100.csv"

# Skip the first line and read the file line by line
tail -n +2 "$combos_file" | while IFS=, read -r maf popmaf distance d r2; do

  for n_snps in $n_snps_per_str; do

    for n_snp_rep in $(seq 1 $n_snp_reps); do

      # Make string for combo and snplist_id
      combo=${maf}_${popmaf}_${distance}_${d}_${r2}
      snplist_id=${experiment}_${combo}_${n_snps}_${n_snp_rep}

      # Submit the record matching pipeline
      for n_ind_rep in $(seq 1 $n_ind_reps); do

        n_run=${snplist_id}_${n_ind_rep}

        if [[ ! -f "$out_dir/run_$n_run/match_accuracies.csv" ]]; then
          submit_script="$submit_dir/submit_${n_run}.sh"
          rm -f "$submit_script"

          cat <<EOT > "$submit_script"
#!/bin/bash
#SBATCH --job-name=rm_${n_run}
#SBATCH --time=$time_for_job
#SBATCH --mem=2G
#SBATCH --output=logs_${fraction}/rm_${n_run}.log
#SBATCH --error=logs_${fraction}/rm_${n_run}.err
#SBATCH -p hns,owners,normal
#SBATCH -c 1

$DIR_CODE/record_matching_pipeline_combinations.sh \
$snplist_id $n_ind_rep $fraction

EOT

          echo "Created and submitting $submit_script for $n_run"
          sbatch "$submit_script" || { echo "Error submitting $submit_script"; exit 1; }
        fi

      done
    done
  done
done
