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

  maf = as.numeric(args[4])
  popmaf = as.numeric(args[5])
  distance = as.numeric(args[6])
  payseur_d = as.numeric(args[7])
  ld_r2 = as.numeric(args[8])

}

combination = str_interp("${maf}_${popmaf}_${distance}_${payseur_d}_${ld_r2}")
snplist_id = str_interp("${experiment}_${combination}_${n_snps}_${n_snp_rep}")

set.seed(n_snp_rep)


#------------------#
# SOURCE VARIABLES #
#------------------#
dir_data_processed = Sys.getenv("DIR_DATA_PROCESSED")
dir_data_snp_lists = Sys.getenv("DIR_DATA_SNP_LISTS")
dir_tmp = Sys.getenv("DIR_TMP")
codis_strs = unlist(str_split(Sys.getenv("CODIS_STRS"), " "))


#------------#
# FUNCTIONS  #
#------------#
filter_maf = function(df, maf=NA){

  if (is.na(maf) == FALSE){
    df = df %>%
      filter(MAF > maf)
  }

return(df)

}

filter_popmaf = function(df, popmaf=NA){

  if (is.na(popmaf) == FALSE){
    df = df %>%
      filter(AFR_MAF > popmaf, AMR_MAF > popmaf,
             EAS_MAF > popmaf, EUR_MAF > popmaf,
             SAS_MAF > popmaf)
  }

  return(df)

}

filter_distance = function(df, dist=NA){

  if (is.na(dist) == FALSE){
    df = df %>%
      filter(abs(dist_to_STR) < dist)
  }

  return(df)

}

filter_payseur_d = function(df, payseur_d=NA){

  if (is.na(payseur_d) == FALSE){
    df = df %>%
      filter(Payseur_d > payseur_d)
  }

  return(df)

}

filter_ld_pruning = function(df, str, r2=NA){

  if (is.na(r2) == FALSE){
  ld_file_path = str_interp("${dir_tmp}/${str}_r2_${r2}.prune.in")
  ld_df = read.csv(ld_file_path, header=FALSE)
  names(ld_df) = 'snp_id'

  df = df %>%
    filter(SNP %in% ld_df$snp_id)

  }

  return(df)

}


#-----------------#
# MAKE SNP LISTS  #
#-----------------#
setwd(dir_data_snp_lists)

for (str in codis_strs){

  snp_file_path = str_interp("${dir_data_processed}/${str}_snps.csv")
  df = read.csv(snp_file_path)

  # Perform combinatory filtering
  df = filter_maf(df, maf)
  df = filter_popmaf(df, popmaf)
  df = filter_distance(df, distance)
  df = filter_payseur_d(df, payseur_d)
  df = filter_ld_pruning(df, str, ld_r2)

  # Test
  test_snps_df = df %>%
    filter(SNP != str) %>%
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
