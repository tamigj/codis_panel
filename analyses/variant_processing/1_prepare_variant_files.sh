#!/bin/bash

codis_str=$1

source /scratch/groups/noahr/tami/codis_panel/config.sh
ml R/4.2

cd $DIR_PLINK

  echo ${codis_str}

  # Calculate MAF
  ./plink \
    --vcf "$DIR_DATA_PROCESSED/${codis_str}_withSTR_GT.vcf" \
    --make-bed \
    --freq \
    --hardy \
    --out "$DIR_TMP/${codis_str}" &&

  # Calculate pop-MAF sequentially within the same function
  for superpop in $SUPERPOPULATIONS; do

    ./plink \
      --vcf "$DIR_DATA_PROCESSED/${codis_str}_withSTR_GT.vcf" \
      --keep "$DIR_DATA_RAW/${superpop}_ids_plink.txt" \
      --freq \
      --out "$DIR_TMP/${codis_str}_${superpop}"

  done &&

  cd $DIR_ANALYSES_VARIANT_PROCESSING &&

  # Calculate Payseur LD
  Rscript 1a_calculate_payseur_ds_for_str.R "${codis_str}" &&

  # Merge all
  Rscript 1b_merge_bim_popMAF_and_LD.R "${codis_str}"
