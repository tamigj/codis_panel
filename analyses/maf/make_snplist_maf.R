#!/usr/bin/env Rscript

library(dplyr)
library(stringr)
library(data.table)

#--------------------#
# PROCESS ARGUMENTS  #
#--------------------#
args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
  stop("You must specify the STR", call.=FALSE)
} else if (length(args) != 0) {
  experiment = as.character(args[1])
  n_snps = as.numeric(args[2])
  n_snp_rep = as.character(args[3])

  maf_filter = as.numeric(args[4])
}

snplist_id=str_interp("${experiment}_${maf_filter}_${n_snps}_${n_snp_rep}")

set.seed(n_snp_rep)


#------------------#
# SOURCE VARIABLES #
#------------------#
dir_data_processed = Sys.getenv("DIR_DATA_PROCESSED")
dir_data_snp_lists = Sys.getenv("DIR_DATA_SNP_LISTS")
codis_strs = unlist(str_split(Sys.getenv("CODIS_STRS"), " "))
n_snps_per_str = as.numeric(unlist(str_split(Sys.getenv("ALL_PANEL_SIZES"), " ")))


#--------------------#
# PREPARE SNP LISTS  #
#--------------------#
setwd(dir_data_snp_lists)

for (str in codis_strs){

  snp_file_path = str_interp("${dir_data_processed}/${str}_snps.csv")
  snp_df = read.csv(snp_file_path)

  # Test
  test_snps_df = snp_df %>%
    filter(SNP != str) %>%
    filter(MAF >= maf_filter) %>%
    sample_n(n_snps) %>%
    select(SNP)

  # Reference
  ref_snps_df = rbind(test_snps_df, str)

  # Save files
  ref_path = str_interp("reference_${snplist_id}_${str}.csv")
  test_path = str_interp("test_${snplist_id}_${str}.csv")

  write.csv(ref_snps_df, ref_path, row.names=FALSE, quote=FALSE)
  write.csv(test_snps_df, test_path, row.names=FALSE, quote=FALSE)

}
