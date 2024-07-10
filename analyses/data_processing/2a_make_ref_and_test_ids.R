#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(stringr)

#--------------------#
# PROCESS ARGUMENTS  #
#--------------------#

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
  stop("You must specify the fraction of individuals that should be in the reference list", call.=FALSE)

} else if (length(args) != 0) {

  fraction_ref = as.numeric(args[1])
  run_number = as.numeric(args[2])
  out_dir = as.character(args[3])

}

set.seed(run_number)


#------------------#
# SOURCE VARIABLES #
#------------------#
# Retrieve the environment variables in R
dir_data_raw = Sys.getenv("DIR_DATA_RAW")
dir_data_processed = Sys.getenv("DIR_DATA_PROCESSED")

#------------#
# LOAD DATA  #
#------------#
setwd(dir_data_raw)

vcf_df = data.frame(fread('CSF1PO_halfwindow500000WithSTR.vcf', skip=10))


#--------------#
# EXTRACT IDs  #
#--------------#
ids_1kpg = names(vcf_df %>% select(starts_with("HG") | starts_with("NA")))

n_indiv = length(ids_1kpg)
n_test_indiv = fraction_ref*n_indiv


#---------------------------------#
# DEFINE REFERENCE AND TEST SETS  #
#---------------------------------#
ref_ids = ids_1kpg[sample(n_indiv, n_test_indiv)]
test_ids = ids_1kpg[!(ids_1kpg %in% ref_ids)]

ref_ids_df = data.frame(ref_ids)
test_ids_df = data.frame(test_ids)

n_ref_indiv = n_indiv-n_test_indiv

print(str_interp("There are ${n_ref_indiv} reference (${fraction_ref*100}%) and ${n_test_indiv} test individuals."))


#-------------------------#
# SAVE INDIVIDUAL LISTS   #
#-------------------------#
setwd(out_dir)

write.csv(ref_ids_df, str_interp("ref_ids_${run_number}.csv"), row.names=FALSE, quote=FALSE)
write.csv(test_ids_df, str_interp("test_ids_${run_number}.csv"), row.names=FALSE, quote=FALSE)

print(str_interp("Saved reference file at ${out_dir}/ref_ids_${run_number}.csv"))
print(str_interp("Saved test file at ${out_dir}/test_ids_${run_number}.csv"))
