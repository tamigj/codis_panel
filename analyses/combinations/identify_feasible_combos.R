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

maf_filters = c(as.numeric(unlist(str_split(Sys.getenv("MAF_FILTERS"), " "))), NA)
popmaf_filters = c(as.numeric(unlist(str_split(Sys.getenv("POPMAF_FILTERS"), " "))), NA)
distance_filters = c(as.numeric(unlist(str_split(Sys.getenv("DISTANCE_FILTERS"), " "))), NA)
payseur_d_filters = c(as.numeric(unlist(str_split(Sys.getenv("PAYSEUR_D_FILTERS"), " "))), NA)


#------------#
# FUNCTIONS  #
#------------#
filter_maf = function(df, maf=NA){

  if (is.na(maf) == FALSE){
    df = df %>%
      filter(MAF >= maf)
  }

  return(df)

}

filter_popmaf = function(df, popmaf=NA){

  if (is.na(popmaf) == FALSE){

    if (as.numeric(popmaf) == 0){
      df = df %>%
        filter(AFR_MAF > popmaf, AMR_MAF > popmaf,
               EAS_MAF > popmaf, EUR_MAF > popmaf,
               SAS_MAF > popmaf)
    } else {
      df = df %>%
        filter(AFR_MAF >= popmaf, AMR_MAF >= popmaf,
               EAS_MAF >= popmaf, EUR_MAF >= popmaf,
               SAS_MAF >= popmaf)
    }

  }

  return(df)

}

filter_distance = function(df, dist=NA){

  if (is.na(dist) == FALSE){
    df = df %>%
      filter(abs(dist_to_STR) <= dist)
  }

  return(df)

}

filter_payseur_d = function(df, payseur_d=NA){

  if (is.na(payseur_d) == FALSE){
    df = df %>%
      filter(Payseur_d >= payseur_d)
  }

  return(df)

}


#-----------------#
# MAKE SNP LISTS  #
#-----------------#
summary_df = data.frame(matrix(nrow=0, ncol=6))
names(summary_df) = c('MAF', 'POP_MAF', 'Distance',
                      'Payseur_d', 'STR', 'N')

for (maf in maf_filters){
  for (popmaf in popmaf_filters){
    for (distance in distance_filters){
      for (payseur_d in payseur_d_filters){

          filters = c(maf, popmaf, distance, payseur_d)
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

                n = nrow(df)

                summary_df[nrow(summary_df)+1,] = c(maf, popmaf, distance,
                                                    payseur_d, str, n)


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
  group_by(MAF, POP_MAF, Distance, Payseur_d) %>%
  mutate(N_total = sum(as.numeric(N)),
         min_N = min(as.numeric(N))) %>%
  select(-STR, -N) %>%
  distinct()

feasible_combos_100_df = summary_total_snps_df %>%
  filter(min_N >= 100) %>%
  select(MAF, POP_MAF, Distance, Payseur_d)


# Save data
write.csv(summary_df,
          str_interp('${dir_output_sumstats}/combinations_n_snps_all_strs.csv'),
          row.names=FALSE, quote=FALSE)

write.csv(summary_total_snps_df,
          str_interp('${dir_output_sumstats}/combinations_n_snps.csv'),
          row.names=FALSE, quote=FALSE)

write.csv(feasible_combos_100_df,
          str_interp('${dir_output_sumstats}/feasible_combinations_100.csv'),
          row.names=FALSE, quote=FALSE)

# Save pop-MAF = 0% combinations specifically
write.csv(feasible_combos_100_df %>% filter(POP_MAF == 0),
          str_interp('${dir_output_sumstats}/feasible_combinations_100_pop_maf.csv'),
          row.names=FALSE, quote=FALSE)
