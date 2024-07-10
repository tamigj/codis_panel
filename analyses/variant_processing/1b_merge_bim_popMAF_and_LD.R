#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(stringr)

#--------------------#
# PROCESS ARGUMENTS  #
#--------------------#
args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
  stop("You must specify the STR", call.=FALSE)
} else if (length(args) != 0) {
  str = as.character(args[1])
}


#------------------#
# SOURCE VARIABLES #
#------------------#
# Retrieve the environment variables in R
dir_tmp = Sys.getenv("DIR_TMP")
dir_data_processed = Sys.getenv("DIR_DATA_PROCESSED")
dir_output_sumstats = Sys.getenv("DIR_OUTPUT_SUMSTATS")


#-----------#
# VARIABLES #
#-----------#
bim_path = str_interp('${dir_tmp}/${str}.bim')
freq_path = str_interp('${dir_tmp}/${str}.frq')
hwe_path = str_interp('${dir_tmp}/${str}.hwe')
ld_path = str_interp('${dir_output_sumstats}/${str}_payseur_ld.csv')

eur_path = str_interp('${dir_tmp}/${str}_EUR.frq')
eas_path = str_interp('${dir_tmp}/${str}_EAS.frq')
sas_path = str_interp('${dir_tmp}/${str}_SAS.frq')
amr_path = str_interp('${dir_tmp}/${str}_AMR.frq')
afr_path = str_interp('${dir_tmp}/${str}_AFR.frq')

out_path = str_interp('${dir_data_processed}/${str}_snps.csv')


#------------#
# LOAD DATA  #
#------------#
bim_df = data.frame(fread(bim_path, header=FALSE))
names(bim_df) = c('CHR', 'SNP', 'RANDOM', 'POS', 'A1', 'A2')

freq_df = data.frame(fread(freq_path, header=TRUE))

hwe_df = data.frame(fread(hwe_path, header=TRUE))
hwe_df = hwe_df %>%
  select(SNP, O.HET.) %>%
  rename(HET = O.HET.)

ld_df = data.frame(fread(ld_path, header=TRUE))

freq_afr_df = data.frame(fread(afr_path, header=TRUE))
freq_amr_df = data.frame(fread(amr_path, header=TRUE))
freq_eas_df = data.frame(fread(eas_path, header=TRUE))
freq_eur_df = data.frame(fread(eur_path, header=TRUE))
freq_sas_df = data.frame(fread(sas_path, header=TRUE))

freq_afr_df = freq_afr_df %>%
  select(SNP, MAF) %>%
  rename(AFR_MAF = MAF)

freq_amr_df = freq_amr_df %>%
  select(SNP, MAF) %>%
  rename(AMR_MAF = MAF)

freq_eas_df = freq_eas_df %>%
  select(SNP, MAF) %>%
  rename(EAS_MAF = MAF)

freq_eur_df = freq_eur_df %>%
  select(SNP, MAF) %>%
  rename(EUR_MAF = MAF)

freq_sas_df = freq_sas_df %>%
  select(SNP, MAF) %>%
  rename(SAS_MAF = MAF)


#------------#
# MERGE DATA #
#------------#
df = merge(bim_df, freq_df, by=c('CHR', 'SNP', 'A1', 'A2'))
df = df %>%
  select(CHR, SNP, A1, A2, POS, MAF)

df = Reduce(function(x, y) left_join(x, y, by = "SNP"),
            list(df, hwe_df, freq_afr_df, freq_amr_df, freq_eas_df, freq_eur_df, freq_sas_df))

df = merge(df, ld_df, by='SNP', all.x=TRUE)


#----------------------#
# ADD DISTANCE TO STR  #
#----------------------#
pos = df[which(df$SNP == str), 'POS']
df$dist_to_STR = df$POS - pos

write.csv(df, out_path, row.names=FALSE, quote=FALSE)

print(str_interp("I saved the file at ${out_path}."))
