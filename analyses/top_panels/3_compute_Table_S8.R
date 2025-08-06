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
top_panels = c('combinations_0.05_NA_125000_NA',
               'combinations_0.1_NA_125000_NA',
               'combinations_NA_0_62500_NA',
               'combinations_NA_0_125000_NA',
               'combinations_NA_0.01_62500_NA',
               'combinations_NA_0.05_125000_NA')

top_panels_txt = c("MAF≥5%, Distance≤0.125Mb",
                   "MAF≥10%, Distance≤0.125Mb",
                   "Pop-MAF≥0%, Distance≤0.0625Mb",
                   "Pop-MAF≥0%, Distance≤0.125Mb",
                   "Pop-MAF≥1%, Distance≤0.0625Mb",
                   "Pop-MAF≥5%, Distance≤0.125Mb")

top_panels_df = data.frame(cbind(top_panels, top_panels_txt))
names(top_panels_df) = c('code', 'text')

panel_sizes = c(50, 100)


#------------#
# FUNCTIONS  #
#------------#
organize_results_in_table = function(results_df, top_panels_df){
  
  results_df$Combination = 
    top_panels_df$text[match(results_df$Combination, top_panels_df$code)]
  
  clean_df = results_df %>%
    group_by(Ratio, `Minimum match score`, Combination, `Panel size`) %>%
    summarize(sum_numbers = sum(as.numeric(Number)),
              n = n(),
              total = n*626,
              Percentage = sum_numbers/total) %>%
    ungroup() %>%
    
    mutate(Percentage = round(Percentage, 2)) %>%
    select(-sum_numbers, -n, -total) %>%
    
    mutate(Merged_panels = paste(Combination, `Panel size`, sep="_"),
           Merged_ratios = paste("x", Ratio, `Minimum match score`, sep='_')) %>%
    
    select(-Combination, -`Panel size`, -Ratio, -`Minimum match score`) %>%
    
    spread(key=Merged_ratios, value=Percentage) %>%
    separate(Merged_panels, into=c("Panel", "Panel size"), sep='_') %>%
    arrange(Panel, desc(`Panel size`))
  
  clean_df[nrow(clean_df)+1, ] = NA
  clean_df[nrow(clean_df)+1, ] = NA
  clean_df = data.frame(clean_df[c(13:14, 1:12), ])
  
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
           10^11, 10^12, 10^13, 10^14, 10^15, 10^16, 10^17)

results_df = data.frame(matrix(nrow=0, ncol=6))
names(results_df) = c("Ratio", "Minimum match score",
                      "Combination", "Panel size", 
                      "Number", "Percent")


for (combination in top_panels){
  for (panel_size in panel_sizes){
    
    panel = str_interp("${combination}_${panel_size}")
    print(panel)
    
    for (i in 1:10){
      for (j in 1:10){
        
        dir = str_interp("/scratch/groups/noahr/tami/codis_panel/output/experiments/0.75/run_${panel}_${i}_${j}")
        
        # Process MSM if file exists 
        if (file.exists(dir) & file.exists(str_interp("${dir}/match_score_matrix.csv"))){
          
          setwd(dir)
          
          df_msm = data.frame(fread("match_score_matrix.csv"))
          df_msm = df_msm %>% select(-V1)
          
          for (ratio in ratios){
            
            min_match_score = log(ratio)  
            
            n_exceeding = sum(diag(as.matrix(df_msm)) > min_match_score)
            perc_exceeding = n_exceeding/nrow(df_msm)
            
            result = c(ratio, round(min_match_score, 2),
                       combination, panel_size*18,
                       n_exceeding, round(perc_exceeding, 2))
            
            results_df[nrow(results_df)+1, ] = result
          }
        }
        
      }
    }
    
  }
}

clean_df = organize_results_in_table(results_df, top_panels_df)

#-------------------------------------------------------------------------------
#------------#
# VARIABLES  #
#------------#
top_panels = c("random_all", "random_15str_all")
top_panels_txt = c("All SNPs (18 STRs)", "All SNPs (15 STRs)")
top_panels_sizes = c(192672, 161968)

top_panels_df = data.frame(cbind(top_panels, top_panels_txt, top_panels_sizes))
names(top_panels_df) = c('code', 'text', 'panel_size')


#-----------#
# ANALYSIS  #
#-----------#
all_results_df = data.frame(matrix(nrow=0, ncol=6))
names(all_results_df) = c("Ratio", "Minimum match score",
                          "Combination", "Panel size", 
                          "Number", "Percent")

for (panel in top_panels){
  
  cat("Processing panel:", panel, "\n")
  
  # Get the correct panel size for this panel
  panel_size = top_panels_df$panel_size[top_panels_df$code == panel]
  
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
          
          result = c(ratio, round(min_match_score, 2),
                     panel, panel_size,
                     n_exceeding, round(perc_exceeding, 2))
          
          all_results_df[nrow(all_results_df)+1, ] = result
        }
      }
    }
  }
}

clean_df_all = organize_results_in_table(all_results_df, top_panels_df)

resulting_df = rbind(clean_df, clean_df_all[4:3,])


#-------------#
# SAVE TABLE  #
#-------------#
setwd("/scratch/groups/noahr/tami/codis_panel/output/rm_summaries/0.75/")
write.csv(resulting_df, "fraction_true_matches_TableS8.csv", row.names=FALSE)
