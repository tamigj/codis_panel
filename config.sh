#!/bin/bash

# CODIS STRs
export CODIS_STRS="CSF1PO D10S1248 D12S391 D13S317 D18S51 D19S433 D1S1656 D22S1045 D2S1338 D2S441 D3S1358 D5S818 D7S820 D8S1179 FGA TH01 TPOX vWA"

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
export DIR_CODE="${DIR_CODIS}/code/"

# Analyses directories
export DIR_ANALYSES="${DIR_CODIS}/analyses/"
export DIR_ANALYSES_VARIANT_PROCESSING="${DIR_CODIS}/analyses/variant_processing"

# Software directories
export DIR_SOFTWARE="${DIR_CODIS}/software"
export DIR_VCFTOOLS="${DIR_CODIS}/software/vcftools_0.1.13/bin/"
export DIR_PLINK="${DIR_CODIS}/software/plink/"

# Tmp directory
export DIR_TMP="${DIR_CODIS}/tmp"

# Output directory
export DIR_OUTPUT="${DIR_CODIS}/output"
export DIR_OUTPUT_SUMSTATS="${DIR_CODIS}/output/sumstats"
export DIR_OUTPUT_EXPERIMENTS="${DIR_CODIS}/output/experiments"
export DIR_OUTPUT_RM_SUMMARIES="${DIR_CODIS}/output/rm_summaries"

# Record matching variables
export CODIS_STRS_RM="CSF1PO,D10S1248,D12S391,D13S317,D18S51,D19S433,D1S1656,D22S1045,D2S1338,D2S441,D3S1358,D5S818,D7S820,D8S1179,FGA,TH01,TPOX,vWA"
export BEAGLE_JAR="${DIR_CODIS}/software/beagle.22Jul22.46e.jar"
export VCF_EXE="${DIR_CODIS}/software/vcftools_0.1.13/bin/vcf"
