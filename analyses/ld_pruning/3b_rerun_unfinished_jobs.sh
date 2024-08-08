#!/bin/bash

#SBATCH --job-name=rerun_jobs_0.5
#SBATCH --time=00:30:00
#SBATCH --mem=1G
#SBATCH --output=logs_0.5/rerun_jobs_0.5.log
#SBATCH --error=logs_0.5/rerun_jobs_0.5.err
#SBATCH -p hns,owners,normal
#SBATCH -c 1

source /scratch/groups/noahr/tami/codis_panel/config.sh

#---------------------------#
# EDIT FOR EACH EXPERIMENT  #
#---------------------------#
experiment='pruning'
n_snps_per_str=$SELECTED_PANEL_SIZES
n_snp_reps=10
n_ind_reps=10
fraction=0.5
filters=$R2_THRESHOLDS

time_for_job="15:00:00"

# Directory paths
submit_dir="submit_$fraction"
out_dir="$DIR_OUTPUT_EXPERIMENTS/$fraction"

# Main loop
for filter in $filters; do
    for n_snps in $n_snps_per_str; do
        for n_snp_rep in $(seq 1 $n_snp_reps); do
            snplist_id=${experiment}_${filter}_${n_snps}_${n_snp_rep}

            for n_ind_rep in $(seq 1 $n_ind_reps); do
                n_run=${snplist_id}_${n_ind_rep}

                if [[ ! -f "$out_dir/run_$n_run/match_accuracies.csv" ]]; then
                    submit_script="$submit_dir/submit_${n_run}.sh"
                    rm -f $submit_script

                    cat <<EOT > "$submit_script"
#!/bin/bash
#SBATCH --job-name=rm_${n_run}
#SBATCH --time=$time_for_job
#SBATCH --mem=2G
#SBATCH --output=logs_${fraction}/rm_${n_run}.log
#SBATCH --error=logs_${fraction}/rm_${n_run}.err
#SBATCH -p hns,owners,normal
#SBATCH -c 1

$DIR_CODE/record_matching_pipeline_variant_characteristics.sh \\
$experiment $n_snps $n_snp_rep $n_ind_rep $fraction $filter
EOT
                    echo "Created and submitting $submit_script for $n_run"
                    sbatch "$submit_script"
                fi
            done
        done
    done
done
