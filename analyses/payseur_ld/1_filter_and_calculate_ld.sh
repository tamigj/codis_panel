#!/bin/bash

#SBATCH --job-name=payseur_snplists_0.75
#SBATCH --time=00:30:00
#SBATCH --mem=1G
#SBATCH --output=logs/payseur_snplists_0.75.log
#SBATCH --error=logs/payseur_snplists_0.75.err
#SBATCH -p hns,owners,normal
#SBATCH -c 1

source /scratch/groups/noahr/tami/codis_panel/config.sh

#---------------------------#
# EDIT FOR EACH EXPERIMENT  #
#---------------------------#
experiment='payseur_d'
n_snps_per_str=$SELECTED_PANEL_SIZES
n_snp_reps=10
fraction=0.75

###############################################################
ml R/4.2

mkdir -p logs_${fraction}
mkdir -p submit_${fraction}

for payseur_d_filter in $PAYSEUR_D_FILTERS;
do

  for n_snps in $n_snps_per_str;
  do

    for n_snp_rep in $(seq 1 $n_snp_reps);
    do

      snplist_id=${experiment}_${payseur_d_filter}_${n_snps}_${n_snp_rep}

      # 1. Prepare SNP list
      Rscript $DIR_ANALYSES_PAYSEUR/make_snplist_payseur_d.R ${experiment} ${n_snps} ${n_snp_rep} ${payseur_d_filter}

      # 2. Calculate LD
      cd $DIR_PLINK

      for str in $CODIS_STRS;
      do
        for r2 in $R2_THRESHOLDS;
        do

          ./plink \
             --vcf $DIR_DATA_PROCESSED/${str}_withSTR_GT.vcf \
             --extract $DIR_DATA_SNP_LISTS/test_payseur_d_${payseur_d_filter}_${n_snps}_${n_snp_rep}_${str}.csv \
             --r2 \
             --make-bed \
             --out $DIR_OUTPUT_PAYSEUR_LD/payseur_d_${payseur_d_filter}_${n_snps}_${n_snp_rep}_${str}

        done
      done

      rm $DIR_OUTPUT_PAYSEUR_LD/*bed
      rm $DIR_OUTPUT_PAYSEUR_LD/*bim
      rm $DIR_OUTPUT_PAYSEUR_LD/*fam
      rm $DIR_OUTPUT_PAYSEUR_LD/*nosex

    done
  done
done
