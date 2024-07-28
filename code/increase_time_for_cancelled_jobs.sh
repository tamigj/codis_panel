# Directory containing the error files
err_dir=$1

# Directory containing the submit scripts
submit_dir=$2

# Time settings
old_time=$3
new_time=$4

# Search for CANCELLED jobs in the error files
grep "CANCELLED" $err_dir/*.err | while read -r line; do
    # Extract the job identifier from the error file name
    job_file=$(echo "$line" | cut -d: -f1)
    job_id=$(basename "$job_file" .err | cut -d_ -f2-)

    # Locate the corresponding submit script
    submit_script="$submit_dir/submit_${job_id}.sh"

    if [[ -f "$submit_script" ]]; then
        # Update the time limit in the submit script
        sed -i "s/#SBATCH --time=$old_time/#SBATCH --time=$new_time/" "$submit_script"
        echo "Updated $submit_script"
        
        # Print the updated submit script for verification
        sbatch $submit_script
    else
        echo "Submit script $submit_script not found"
    fi
done
