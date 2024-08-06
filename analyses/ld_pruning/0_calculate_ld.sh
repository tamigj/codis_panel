#!/bin/bash

#SBATCH --job-name=logs/calc_ld_r2
#SBATCH --time=01:00:00
#SBATCH --mem=1G
#SBATCH --output=logs/calc_ld_r2.log
#SBATCH --error=logs/calc_ld_r2.err
#SBATCH -p normal,owners,hns
#SBATCH -c 1

mkdir -p logs
source /scratch/groups/noahr/tami/codis_panel/config.sh

cd $DIR_PLINK

for str in $CODIS_STRS;
do

  echo ${codis_snp} >$DIR_TMP/exclude_${str}.csv

  ./plink \
     --vcf $DIR_DATA_PROCESSED/${str}_withSTR_GT.vcf \
     --exclude $DIR_TMP/exclude_${str}.csv \
     --make-bed \
     --out $DIR_TMP/${str}

  for r2 in $R2_THRESHOLDS;
  do

     ./plink \
       --bfile $DIR_TMP/${str} \
       --indep-pairwise 100 10 ${r2} \
       --out $DIR_TMP/${str}_r2_${r2}
  done

  rm $DIR_TMP/${str}.bed
  rm $DIR_TMP/${str}.bim
  rm $DIR_TMP/${str}.fam

done
