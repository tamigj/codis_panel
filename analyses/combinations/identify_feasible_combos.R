#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)


#------------#
# VARIABLES  #
#------------#
codis_strs = unlist(str_split(Sys.getenv("CODIS_STRS"), " "))
dir_data_processed = Sys.getenv("DIR_DATA_PROCESSED")
dir_tmp = Sys.getenv("DIR_TMP")
dir_output_sumstats = Sys.getenv("DIR_OUTPUT_SUMSTATS")

maf_filters = c(0.01, NA)
popmaf_filters = c(0.01, NA)
distance_filters = c(250000, NA)
payseur_d_filters = c(0.3, NA)
r2_filters = c(0.7, NA)


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
summary_df = data.frame(matrix(nrow=0, ncol=7))
names(summary_df) = c('MAF', 'POP_MAF', 'Distance',
                      'Payseur_d', 'R2', 'STR', 'N')

for (maf in maf_filters){
  for (popmaf in popmaf_filters){
    for (distance in distance_filters){
      for (payseur_d in payseur_d_filters){
        for (r2 in r2_filters){

          filters = c(maf, popmaf, distance, payseur_d, r2)
          n_filters = length(filters) - sum(is.na(filters))
          n_maf_filters = 2 - sum(is.na(c(maf, popmaf)))

          # Combination of filters must be >2 & must not have both MAF and pop-MAF
          if (n_filters > 1){
            if (n_maf_filters < 2){

              # Perform all filtering
              for (str in codis_strs){

                snp_file_path = str_interp("${dir_data_processed}/${str}_snps.csv")
                df = read.csv(snp_file_path)

                df = filter_maf(df, maf)
                df = filter_popmaf(df, popmaf)
                df = filter_distance(df, distance)
                df = filter_payseur_d(df, payseur_d)
                df = filter_ld_pruning(df, str, r2)

                n = nrow(df)

                summary_df[nrow(summary_df)+1,] = c(maf, popmaf, distance,
                                                    payseur_d, r2, str, n)

              }

            }
          }

        }
      }
    }
  }
}

# Summarize data
summary_total_snps_df = summary_df %>%
  distinct() %>%
  group_by(MAF, POP_MAF, Distance, Payseur_d, R2) %>%
  mutate(N_total = sum(as.numeric(N)),
         min_N = min(as.numeric(N))) %>%
  select(-STR, -N) %>%
  distinct()

# Feasible combinations
feasible_combos_df = summary_total_snps_df %>%
  filter(min_N >= 100) %>%
  select(MAF, POP_MAF, Distance, Payseur_d, R2)


# Save data
write.csv(summary_df,
          str_interp('${dir_output_sumstats}/combinations_n_snps_all_strs.csv'),
          row.names=FALSE, quote=FALSE)

write.csv(summary_total_snps_df,
          str_interp('${dir_output_sumstats}/combinations_n_snps.csv'),
          row.names=FALSE, quote=FALSE)

write.csv(feasible_combos_df,
          str_interp('${dir_output_sumstats}/feasible_combinations.csv'),
          row.names=FALSE, quote=FALSE)
