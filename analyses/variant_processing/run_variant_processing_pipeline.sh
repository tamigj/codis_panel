#!/bin/bash

#SBATCH --job-name=submit_variant_processing_pipeline
#SBATCH --time=00:30:00
#SBATCH --mem=1G
#SBATCH --output=submit_variant_processing_pipeline.log
#SBATCH --error=submit_variant_processing_pipeline.err
#SBATCH -p hns,owners
#SBATCH -c 1

# Source the configuration file
source /scratch/groups/noahr/tami/codis_panel/config.sh

# Loop through each codis_str
for codis_str in $CODIS_STRS; do

  # Create the submit script
  cat <<EOT > submit_${codis_str}.sh
#!/bin/bash
#SBATCH --job-name=variant_processing_${codis_str}
#SBATCH --time=15:00:00
#SBATCH --mem=10G
#SBATCH --output=variant_processing_${codis_str}.log
#SBATCH --error=variant_processing_${codis_str}.err
#SBATCH -p hns,owners
#SBATCH -c 1

./1_prepare_variant_files.sh ${codis_str}
EOT

  # Make the script executable
  chmod +x submit_${codis_str}.sh

  # Submit the job script
  sbatch submit_${codis_str}.sh

done
