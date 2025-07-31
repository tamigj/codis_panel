library(stringr)
library(dplyr)
library(data.table)
library(RecordMatching)
library(tidyr)


#------------#
# VARIABLES  #
#------------#
top_panels = c("random_15str_all")
top_panels_txt = ("All SNPs (15 STRs)")


#-----------#
# ANALYSIS  #
#-----------#
results_df = data.frame(matrix(nrow=0, ncol=6))
names(results_df) = c("Ratio", "Minimum match score",
                      "Combination", "Panel size", 
                      "Number", "Percent")

df_all = data.frame(matrix(nrow=0, ncol=2))
names(df_all) = c('Diagonal_values', 'Replicate')

for (panel in top_panels){
  
  for (i in 1){
    for (j in 1:100){
      
      dir = str_interp("/scratch/groups/noahr/tami/codis_panel/output/experiments/0.75/run_${panel}_${i}_${j}")
      
      # Process MSM if file exists 
      if (file.exists(dir) & file.exists(str_interp("${dir}/match_score_matrix.csv"))){
        
        setwd(dir)
        
        for (ratio in ratios){
          
          df_msm = data.frame(fread("match_score_matrix.csv"))
          df_msm = df_msm %>% select(-V1)
          
          df_msm_diag = cbind(diag(as.matrix(df_msm)), rep(j, 626))
          
          df_all = rbind(df_all, df_msm_diag)
          
        }
      }
      
    }
  }
  
}

names(df_all) = c("Diagonal_values", 'Replicate')

ggplot(df_all, aes(x=Diagonal_values, fill=Replicate)) +
  geom_density(alpha=0.1)
