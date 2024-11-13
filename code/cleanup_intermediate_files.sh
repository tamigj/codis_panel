#!/bin/bash

snplist_id=$1
n_run=$2
fraction=$3

# Directories (data and output)
DIR_DATA_NRUN="$DIR_DATA_EXPERIMENTS/$fraction/run_$n_run"
DIR_OUTPUT_NRUN="$DIR_OUTPUT_EXPERIMENTS/$fraction/run_$n_run"

DIR_DATA_NRUN_REF="$DIR_DATA_NRUN/reference"
DIR_DATA_NRUN_SNP="$DIR_DATA_NRUN/SNP"
DIR_DATA_NRUN_STR="$DIR_DATA_NRUN/STR"

# Delete data files
rm $DIR_DATA_NRUN_REF/*
rm $DIR_DATA_NRUN_SNP/*
rm $DIR_DATA_NRUN_STR/*

# Print confirmation after data deletion
echo "Deleted all reference, SNP and STR plink files in $DIR_DATA_NRUN/"

# Delete intermediate output files
rm $DIR_OUTPUT_NRUN/*FORMAT
rm $DIR_OUTPUT_NRUN/*vcf
rm $DIR_OUTPUT_NRUN/ref_alfrq/*

# Print confirmation after intermediate output deletion
echo "Deleted all plink files in $DIR_OUTPUT_NRUN"
