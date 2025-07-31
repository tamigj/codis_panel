rm(list=ls())

#-----------------#
# LOAD LIBRARIES  #
#-----------------#
library(stringr)
library(dplyr)
library(data.table)
library(RecordMatching)
library(tidyr)


#------------#
# VARIABLES  #
#------------#
top_panels = c("random_all")
top_panels_txt = ("All SNPs (N=192,672)")

top_panels_df = data.frame(cbind(top_panels, top_panels_txt))
names(top_panels_df) = c('code', 'text')


#------------#
# FUNCTIONS  #
#------------#
organize_results_in_table_all_data = function(results_df){
  
  clean_df = results_df %>%
    group_by(Ratio, `Minimum match score`, Combination, `Panel size`) %>%
    summarize(sum_numbers = sum(as.numeric(Number)),
              n = n(),
              total = n*626,
              Percentage = sum_numbers/total) %>%
    ungroup() %>%
    
    mutate(Percentage = round(Percentage, 3)) %>%
    select(-sum_numbers, -n, -total) %>%
    
    mutate(Merged_panels = paste(Combination, `Panel size`, sep="_"),
           Merged_ratios = paste("x", Ratio, `Minimum match score`, sep='_')) %>%
    
    select(-Combination, -`Panel size`, -Ratio, -`Minimum match score`) %>%
    
    spread(key=Merged_ratios, value=Percentage) %>%
    separate(Merged_panels, into=c("Panel", "Panel size"), sep='_') %>%
    arrange(Panel, desc(`Panel size`))
  
  clean_df[nrow(clean_df)+1, ] = NA
  clean_df[nrow(clean_df)+1, ] = NA
  clean_df = data.frame(clean_df[c(2:3, 1), ])
  
  for (i in 3:ncol(clean_df)){
    
    colname = names(clean_df)[i]
    
    ratio = unlist(strsplit(colname, "_"))[2]
    score = unlist(strsplit(colname, "_"))[3]
    
    clean_df[1,i] = ratio
    clean_df[2,i] = score
    
  }
  
  return(clean_df)
  
}


#-----------#
# ANALYSIS  #
#-----------#
ratios = c(10^0, 10^1, 10^2, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9, 10^10,
           10^11, 10^12, 10^13, 10^14, 10^15, 10^16, 10^17, 10^18, 10^19, 10^20)

results_df = data.frame(matrix(nrow=0, ncol=6))
names(results_df) = c("Ratio", "Minimum match score",
                      "Combination", "Panel size", 
                      "Number", "Percent")

  
for (panel in top_panels){

  for (i in 1){
    for (j in 1:100){
          
      dir = str_interp("/scratch/groups/noahr/tami/codis_panel/output/experiments/0.75/run_${panel}_${i}_${j}")
          
      # Process MSM if file exists 
      if (file.exists(dir) & file.exists(str_interp("${dir}/match_score_matrix.csv"))){
            
        setwd(dir)
            
        for (ratio in ratios){
            
          min_match_score = log(ratio)
            
          df_msm = data.frame(fread("match_score_matrix.csv"))
          df_msm = df_msm %>% select(-V1)
              
          n_exceeding = sum(diag(as.matrix(df_msm)) > min_match_score)
          perc_exceeding = n_exceeding/nrow(df_msm)
              
          result = c(ratio, round(min_match_score, 3),
                     panel, 'all data',
                     n_exceeding, round(perc_exceeding, 3))
              
          results_df[nrow(results_df)+1, ] = result
          
        }
      }
      
    }
  }

}


#------------#
# SAVE DATA  #
#------------#
setwd("/scratch/groups/noahr/tami/codis_panel/output/rm_summaries/0.75/")
write.csv(results_df, "ratios_and_percentages.csv", row.names=FALSE)

clean_df = organize_results_in_table_all_data(results_df) 
write.csv(clean_df, "ratios_and_percentages_TableS8B.csv", row.names=FALSE)
