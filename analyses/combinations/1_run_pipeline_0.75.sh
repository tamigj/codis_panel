#!/bin/bash

#SBATCH --job-name=submit_jobs_0.75
#SBATCH --time=08:00:00
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
n_snps_per_str=$ALL_PANEL_SIZES_UP_TO_100
n_snp_reps=10
n_ind_reps=10
fraction=0.75

start_line=2
end_line=3
###############################################################

ml R/4.2

mkdir -p logs_${fraction}
mkdir -p submit_${fraction}

combos_file="$DIR_OUTPUT_SUMSTATS/feasible_combinations_100_pop_maf.csv"

if [ ! -f "$combos_file" ]; then
    Rscript identify_feasible_combos.R
fi

# Read the file line by line
while IFS=, read -r maf popmaf distance d; do

    # Increment line counter
    line_counter=$((line_counter + 1))

    # Process only the first three lines
    if [ "$line_counter" -ge $start_line ] && [ "$line_counter" -le $end_line ]; then

    for n_snps in $n_snps_per_str;
    do

      for n_snp_rep in $(seq 1 $n_snp_reps);
      do

        # Make string for combo and snplist_id
        combo=${maf}_${popmaf}_${distance}_${d}
        snplist_id=${experiment}_${combo}_${n_snps}_${n_snp_rep}

        # 1. Prepare SNP list
        Rscript make_snplist_combos.R ${experiment} ${n_snps} ${n_snp_rep} \
                                      ${maf} ${popmaf} ${distance} ${d}

        # 2. Submit the record matching pipeline
        for n_ind_rep in $(seq 1 $n_ind_reps);
        do

          n_run=${snplist_id}_${n_ind_rep}

          # Submit the record_matching_pipeline
          sbatch_script="submit_${n_run}.sh"

          # Check if the sbatch script already exists
          # if [ ! -f submit_$fraction/$sbatch_script ]; then

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

              # Submit the sbatch script
              sbatch $sbatch_script
            fi
          #fi

        done
      done

    done && mv submit*sh submit_$fraction/
  fi

done < "$combos_file"
