#!/bin/bash

source /scratch/groups/noahr/tami/codis_panel/config.sh

ml R/4.2
Rscript 1a_identify_non_CODIS_STRs.R

cd $DIR_VCFTOOLS

for codis_str in $CODIS_STRS;
do

  ./vcftools \
    --vcf $DIR_DATA_RAW/${codis_str}_halfwindow500000WithSTR.vcf \
    --exclude $DIR_TMP/${codis_str}_STRs_to_remove.csv \
    --recode \
    --out $DIR_DATA_PROCESSED/${codis_str}_halfwindow500000WithSTR_bia_only

  awk -v FS="\t" -v OFS="\t" \
    '{for(i=9;i<=NF;i++) {split($i, gt, ":"); $i=gt[1]} print}' \
    $DIR_DATA_PROCESSED/${codis_str}_halfwindow500000WithSTR_bia_only.recode.vcf > \
    $DIR_DATA_PROCESSED/${codis_str}_withSTR_GT.vcf

    echo Done with ${codis_str}.

done &&

rm -r $DIR_TMP/*_STRs_to_remove.csv
rm $DIR_DATA_PROCESSED/*_halfwindow500000WithSTR_bia_only.*
