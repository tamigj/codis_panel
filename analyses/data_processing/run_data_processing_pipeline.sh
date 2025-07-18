#!/bin/bash

#SBATCH --job-name=data_processing_pipeline
#SBATCH --time=10:00:00
#SBATCH --mem=30G
#SBATCH --output=data_processing_pipeline.log
#SBATCH --error=data_processing_pipeline.err
#SBATCH -p hns,owners
#SBATCH -c 1

./0_process_saini_files.sh &&
./1_process_vcfs.sh &&
./2_prepare_partitions.sh &&

mv data_processing_pipeline* logs
