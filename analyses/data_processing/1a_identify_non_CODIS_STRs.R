#!/usr/bin/env Rscript

library(data.table)
library(stringr)


#------------------#
# SOURCE VARIABLES #
#------------------#

# Retrieve the environment variables in R
dir_data_raw = Sys.getenv("DIR_DATA_RAW")
dir_tmp = Sys.getenv("DIR_TMP")
dir_output = Sys.getenv("DIR_OUTPUT")
dir_output_sumstats = Sys.getenv("DIR_OUTPUT_SUMSTATS")
codis_strs = strsplit(Sys.getenv("CODIS_STRS"), " ")[[1]]


#-----------#
# ANALYSIS  #
#-----------#
res_df = data.frame(matrix(nrow=0, ncol=4))
names(res_df) = c('STR', 'n_before', 'n_STRs', 'n_after')

for (codis_str in codis_strs){

  df = fread(str_interp('${dir_data_raw}/${codis_str}_halfwindow500000WithSTR.vcf'))
  tmp_df = df[,3:5]

  for (i in 1:nrow(tmp_df)){
    tmp_df[i,'len_REF'] = nchar(tmp_df[i,'REF'])
    tmp_df[i,'len_ALT'] = nchar(tmp_df[i,'ALT'])
  }

  tmp_MA = tmp_df[which(tmp_df$len_REF != 1 | tmp_df$len_ALT != 1), ]
  tmp_MA = tmp_MA[!(tmp_MA$ID %in% codis_str), ]

  str_list_df = data.frame(tmp_MA$ID)
  names(str_list_df) = c('id')

  ## Record number of SNPs
  n_before = nrow(df)
  n_strs = nrow(str_list_df)
  n_after = n_before-n_strs

  print(str_interp(
    "There are ${n_strs} multi-allelic variants around ${codis_str}."))

  res_df[nrow(res_df)+1, ] = c(codis_str, n_before, n_strs, n_after)

  ## Export file of STRs
  write.csv(str_list_df, str_interp("${dir_tmp}/${codis_str}_STRs_to_remove.csv"),
            row.names=FALSE, quote=FALSE)

}

# Save a file for reference
write.csv(res_df, str_interp("${dir_output_sumstats}/n_STRs_around_CODIS_STRs.csv"),
          row.names=FALSE)
