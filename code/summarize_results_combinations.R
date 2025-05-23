rm(list=ls())

library(data.table)
library(stringr)
library(dplyr)


#--------------------#
# PROCESS ARGUMENTS  #
#--------------------#
args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
  stop("You must specify the fraction of individuals that should be in the reference list", call.=FALSE)

} else if (length(args) != 0) {

  codis_strs = as.character(args[1])      # comma-separated
  experiment = as.character(args[2])
  combos_file = as.character(args[3])
  n_snps_list = as.character(args[4])     # comma-separated

  n_snp_reps = as.numeric(args[5])
  n_ind_reps = as.numeric(args[6])
  fraction = as.numeric(args[7])

  dir_snp_lists = as.character(args[8])
  dir_output_experiments = as.character(args[9])
  dir_output_rm_summaries = as.character(args[10])

}

n_strs = str_count(codis_strs,",") + 1
codis_strs = unlist(strsplit(codis_strs, ",", n_strs))

n_commas = str_count(n_snps_list,",") + 1
n_snps_list = unlist(strsplit(n_snps_list, ",", n_commas))

combos_file_df = read.csv(combos_file)


#----------------------------#
# RECORD MATCHING ACCURACIES #
#----------------------------#
main_df = data.frame(matrix(nrow=0, ncol=10))
names(main_df) = c("combination",
                  "one_to_one", "SNPquery", "STRquery", "needle_in_haystack",
                  "n_snps_per_str", "n_total_snps",
                  "snp_replicate", "ind_replicate")

for (n_snps in n_snps_list){
  for (n_snp_rep in 1:n_snp_reps){
    for (i in 1:nrow(combos_file_df)){

      maf =  combos_file_df[i,'MAF']
      popmaf =  combos_file_df[i,'POP_MAF']
      distance =  combos_file_df[i,'Distance']
      d =  combos_file_df[i,'Payseur_d']

      combo = str_interp("${maf}_${popmaf}_${distance}_${d}")

      # Count total number of SNPs ---------------------------
      snplist_id = str_interp("${experiment}_${combo}_${n_snps}_${n_snp_rep}")
      n_total_snps = 0

      for (str in codis_strs){
        path=str_interp("${dir_snp_lists}/test_${snplist_id}_${str}.csv")
        df_snps = fread(path)
        n = nrow(df_snps)

        n_total_snps = n_total_snps + n
      }

      # Record matching results -------------------------------
      for (n_ind_rep in 1:n_ind_reps){

        n_run = str_interp("${snplist_id}_${n_ind_rep}")

        dir_output_nrun = str_interp("${dir_output_experiments}/${fraction}/run_${n_run}")
        path = str_interp("${dir_output_nrun}/match_accuracies.csv")

        if (file.exists(path)){
          df_experiment = read.csv(path)

          df_experiment$combination = combo
          df_experiment$n_snps_per_str = n_snps
          df_experiment$n_total_snps = n_total_snps
          df_experiment$snp_replicate = n_snp_rep
          df_experiment$ind_replicate = n_ind_rep

          main_df = rbind(main_df, df_experiment)

        } else {
          print(str_interp("No file found for ${n_run}"))
          main_df[nrow(main_df)+1, ] = c(NA, NA, NA, NA,
                                         n_snps, n_total_snps,
                                         n_snp_rep, n_ind_rep)
        }
      }

    }
  }
}


#---------------------#
# SUMMARIZE RESULTS   #
#---------------------#
summary_df = main_df %>%

  group_by(combination, n_snps_per_str) %>%

  mutate(n_replicates = n(),
         experiment = experiment,
         fraction = fraction) %>%

  mutate(one_to_one_median = median(na.omit(as.numeric(one_to_one))),
         one_to_one_min = min(na.omit(as.numeric(one_to_one))),
         one_to_one_max = max(na.omit(as.numeric(one_to_one))),

         SNPquery_median = median(na.omit(as.numeric(SNPquery))),
         SNPquery_min = min(na.omit(as.numeric(SNPquery))),
         SNPquery_max = max(na.omit(as.numeric(SNPquery))),

         STRquery_median = median(na.omit(as.numeric(STRquery))),
         STRquery_min = min(na.omit(as.numeric(STRquery))),
         STRquery_max = max(na.omit(as.numeric(STRquery))),

         needle_in_haystack_median = median(na.omit(as.numeric(needle_in_haystack))),
         needle_in_haystack_min = min(na.omit(as.numeric(needle_in_haystack))),
         needle_in_haystack_max = max(na.omit(as.numeric(needle_in_haystack)))) %>%

  select(experiment, combination, fraction,
         n_snps_per_str, n_total_snps, n_replicates,
         one_to_one_median:needle_in_haystack_max) %>%

  distinct()


#--------------#
# SAVE FILES   #
#--------------#
setwd(str_interp("${dir_output_rm_summaries}/${fraction}"))

write.csv(main_df, str_interp("${experiment}_rm.csv"), row.names=FALSE)
write.csv(summary_df, str_interp("${experiment}_rm_summary.csv"), row.names=FALSE)

print(str_interp("Results saved at ${experiment}_rm.csv"))
