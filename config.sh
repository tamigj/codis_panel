#!/bin/bash

# CODIS STRs
export CODIS_STRS="CSF1PO D10S1248 D12S391 D13S317 D18S51 D19S433 D1S1656 D22S1045 D2S1338 D2S441 D3S1358 D5S818 D7S820 D8S1179 FGA TH01 TPOX vWA"

export CODIS_15_STRS="CSF1PO D13S317 D18S51 D3S1358 D5S818 D7S820 D8S1179 FGA TH01 TPOX vWA D2S441 D10S1248 D19S433 D22S1045"

# Panel sizes (number of SNPs per STR)
export ALL_PANEL_SIZES="10 25 50 75 100 125 250 500 1000 5000"
export ALL_PANEL_SIZES_CSV="10,25,50,75,100,125,250,500,1000,5000"

export SELECTED_PANEL_SIZES="25 50 75 100"
export SELECTED_PANEL_SIZES_CSV="25,50,75,100"

export ALL_PANEL_SIZES_UP_TO_100="10 20 30 40 50 60 70 80 90 100"
export ALL_PANEL_SIZES_UP_TO_100_CSV="10,20,30,40,50,60,70,80,90,100"

# Superpopulations
export SUPERPOPULATIONS="AFR AMR EAS EUR SAS"

# Directories
export DIR_CODIS="/scratch/groups/noahr/tami/codis_panel"

# Data directories
export DIR_DATA="${DIR_CODIS}/data"
export DIR_DATA_RAW="${DIR_CODIS}/data/raw"
export DIR_DATA_PROCESSED="${DIR_CODIS}/data/processed"

export DIR_DATA_PARTITIONS="${DIR_CODIS}/data/partitions"
export DIR_DATA_EXPERIMENTS="${DIR_CODIS}/data/experiments"
export DIR_DATA_SNP_LISTS="${DIR_CODIS}/data/snp_lists"

export DIR_PLINK_MAPS="${DIR_CODIS}/plink.GRCh37.map/"

# Code directory
export DIR_CODE="${DIR_CODIS}/code"
export DIR_CODE_BASELINE="${DIR_CODIS}/code/baseline"
export DIR_CODE_VARIANT_CHARACTERISTICS="${DIR_CODIS}/code/variant_characterstics"

# Analyses directories
export DIR_ANALYSES="${DIR_CODIS}/analyses/"
export DIR_ANALYSES_VARIANT_PROCESSING="${DIR_CODIS}/analyses/variant_processing"
export DIR_ANALYSES_PAYSEUR="${DIR_CODIS}/analyses/payseur_d"

# Software directories
export DIR_SOFTWARE="${DIR_CODIS}/software"
export DIR_VCFTOOLS="${DIR_CODIS}/software/vcftools_0.1.13/bin"
export DIR_PLINK="${DIR_CODIS}/software/plink"

# Tmp directory
export DIR_TMP="${DIR_CODIS}/tmp"

# Output directory
export DIR_OUTPUT="${DIR_CODIS}/output"
export DIR_OUTPUT_SUMSTATS="${DIR_CODIS}/output/sumstats"
export DIR_OUTPUT_EXPERIMENTS="${DIR_CODIS}/output/experiments"
export DIR_OUTPUT_RM_SUMMARIES="${DIR_CODIS}/output/rm_summaries"
export DIR_OUTPUT_PAYSEUR_LD="${DIR_CODIS}/output/payseur_ld"

# Record matching variables
export CODIS_STRS_RM="CSF1PO,D10S1248,D12S391,D13S317,D18S51,D19S433,D1S1656,D22S1045,D2S1338,D2S441,D3S1358,D5S818,D7S820,D8S1179,FGA,TH01,TPOX,vWA"
export CODIS_15_STRS_RM="CSF1PO,D13S317,D18S51,D3S1358,D5S818,D7S820,D8S1179,FGA,TH01,TPOX,vWA,D2S441,D10S1248,D19S433,D22S1045"

export BEAGLE_JAR="${DIR_CODIS}/software/beagle.22Jul22.46e.jar"
export VCF_EXE="${DIR_CODIS}/software/vcftools_0.1.13/bin/vcftools"

# Variant filtering
export MAF_FILTERS="0.01 0.05 0.1"
export POPMAF_FILTERS="0 0.01 0.05"
export DISTANCE_FILTERS="62500 125000 250000"
export PAYSEUR_D_FILTERS="0.3 0.5 0.7"

export MAF_FILTERS_CSV="0.01,0.05,0.1"
export POPMAF_FILTERS_CSV="0,0.01,0.05"
export DISTANCE_FILTERS_CSV="62500,125000,250000"
export PAYSEUR_D_FILTERS_CSV="0.3,0.5,0.7"
