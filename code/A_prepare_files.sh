#!/bin/bash

n_ind_rep=$1
snplist_id=$2
DIR_DATA_NRUN=$3
DIR_DATA_PARTITIONS_FRACTIONS=$4

# Update directories

source /scratch/groups/noahr/tami/codis_panel/config.sh

rm -r $DIR_DATA_NRUN

mkdir $DIR_DATA_NRUN
mkdir "$DIR_DATA_NRUN/reference"
mkdir "$DIR_DATA_NRUN/SNP"
mkdir "$DIR_DATA_NRUN/STR"

for str in $CODIS_STRS;
do

  # Reference ------------------------------------------------------
  "$DIR_VCFTOOLS/vcftools" \
    --vcf "$DIR_DATA_PROCESSED/${str}_withSTR_GT.vcf" \
    --keep "$DIR_DATA_PARTITIONS_FRACTIONS/ref_ids_${n_ind_rep}.csv" \
    --snps "$DIR_DATA_SNP_LISTS/reference_${snplist_id}_${str}.csv" \
    --recode \
    --out "$DIR_DATA_NRUN/reference/${str}"

  # Test (STRs) ----------------------------------------------------
  "$DIR_VCFTOOLS/vcftools" \
    --vcf "$DIR_DATA_PROCESSED/${str}_withSTR_GT.vcf" \
    --keep "$DIR_DATA_PARTITIONS_FRACTIONS/test_ids_${n_ind_rep}.csv" \
    --snp "$str" \
    --extract-FORMAT-info GT \
    --out "$DIR_DATA_NRUN/STR/${str}_STR"

  # Test (SNP) -----------------------------------------------------
  "$DIR_VCFTOOLS/vcftools" \
    --vcf "$DIR_DATA_PROCESSED/${str}_withSTR_GT.vcf" \
    --keep "$DIR_DATA_PARTITIONS_FRACTIONS/test_ids_${n_ind_rep}.csv" \
    --snps "$DIR_DATA_SNP_LISTS/test_${snplist_id}_${str}.csv" \
    --recode \
    --out "$DIR_DATA_NRUN/SNP/${str}_SNP"

done
