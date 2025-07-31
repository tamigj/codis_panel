#!/usr/bin/env Rscript

library(dplyr)
library(stringr)


#--------------------#
# PROCESS ARGUMENTS  #
#--------------------#
args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
  stop("You must specify the STR", call.=FALSE)
} else if (length(args) != 0) {
  snplist_id = as.character(args[1])
}


#------------------#
# SOURCE VARIABLES #
#------------------#
# Retrieve the environment variables in R
dir_data_processed = Sys.getenv("DIR_DATA_PROCESSED")
dir_data_snp_lists = Sys.getenv("DIR_DATA_SNP_LISTS")
codis_strs = unlist(str_split(Sys.getenv("CODIS_15_STRS"), " "))

setwd(dir_data_snp_lists)


#--------------------#
# PREPARE SNP LISTS  #
#--------------------#

for (str in codis_strs){

  snp_file_path = str_interp("${dir_data_processed}/${str}_snps.csv")
  snp_df = read.csv(snp_file_path)

  # Reference
  ref_snps_df = snp_df %>%
    select(SNP)

  # Test
  test_snps_df = snp_df %>%
    select(SNP) %>%
    filter(SNP != str)

  # Save files
  ref_path = str_interp("reference_${snplist_id}_${str}.csv")
  test_path = str_interp("test_${snplist_id}_${str}.csv")

  write.csv(ref_snps_df, ref_path, row.names=FALSE, quote=FALSE)
  write.csv(test_snps_df, test_path, row.names=FALSE, quote=FALSE)

}
