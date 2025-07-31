rm(list=ls())

#-----------------#
# LOAD LIBRARIES  #
#-----------------#
library(stringr)
library(dplyr)
library(data.table)
library(RecordMatching)


#------------#
# VARIABLES  #
#------------#
top_panels = c("random_all")
top_panels_txt = ("All data")

top_panels_df = data.frame(cbind(top_panels, top_panels_txt))
names(top_panels_df) = c('code', 'text')


#------------#
# FUNCTIONS  #
#------------#
clean_up_table = function(results_df){
  
  results_df$Combination = 
    top_panels_df$text[match(results_df$Combination, top_panels_df$code)]
  
  df = results_df %>%
    mutate(Number = as.numeric(Number),
           Percent = as.numeric(Percent)) %>%
    group_by(Ratio, `Minimum match score`, Combination, `Panel size`) %>% 
    summarize(n = n(),
              sum_numbers = sum(Number)) %>%
    mutate(perc = round(sum_numbers/(n*626), 2)) %>%
    select(-n, -sum_numbers) %>%
    ungroup()
  
  # Reorganize table for manuscript
  # Get unique values
  unique_ratios <- unique(df$Ratio)
  unique_min_scores <- unique(df$`Minimum match score`)
  unique_combinations <- unique(df$Combination)
  unique_panel_sizes <- sort(unique((df$`Panel size`)))
  
  # Create all combination-panel pairs, grouped by combination
  combo_panel_pairs <- data.frame()
  for (combo in unique_combinations) {
    for (panel in unique_panel_sizes) {
      combo_panel_pairs <- rbind(combo_panel_pairs, 
                                 data.frame(Combination = combo, 
                                            Panel_size = panel,
                                            stringsAsFactors = FALSE))
    }
  }
  
  # Calculate number of columns needed for perc data
  n_perc_cols <- ceiling(nrow(df) / 12)
  
  # Create the result matrix
  result_matrix <- matrix(NA, nrow = 2 + nrow(combo_panel_pairs), ncol = 2 + n_perc_cols)
  
  # Fill in the first two rows (Ratio and Min match score)
  result_matrix[1, 3:ncol(result_matrix)] <- unique_ratios
  result_matrix[2, 3:ncol(result_matrix)] <- unique_min_scores
  
  # Fill in the combination and panel size columns (rows 3 onwards)
  result_matrix[3:nrow(result_matrix), 1] <- combo_panel_pairs$Combination
  result_matrix[3:nrow(result_matrix), 2] <- combo_panel_pairs$Panel_size
  
  # Fill in the perc values (every 12 elements goes to the next column)
  perc_values <- df$perc
  
  for (i in 1:n_perc_cols) {
    start_idx <- (i - 1) * 12 + 1
    end_idx <- min(i * 12, length(perc_values))
    
    if (start_idx <= length(perc_values)) {
      values_to_fill <- perc_values[start_idx:end_idx]
      # Fill starting from row 3 (since first 2 rows are for ratio and min score)
      result_matrix[3:(2 + length(values_to_fill)), 2 + i] <- values_to_fill
    }
  }
  
  # Convert to data frame
  result_df <- as.data.frame(result_matrix, stringsAsFactors = FALSE)
  
  # Set column names
  colnames(result_df) <- c("Combination", "Panel_size", paste0("Col_", 1:n_perc_cols))
  
  return(result_df)
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

for (ratio in ratios){
  
  min_match_score = log(ratio)
  
  for (panel in top_panels){

      print(panel)
      
      for (i in 1){
        for (j in 1:100){
          
          dir = str_interp("/scratch/groups/noahr/tami/codis_panel/output/experiments/0.75/run_${panel}_${i}_${j}")
          
          # Process MSM if file exists 
          if (file.exists(dir) & file.exists(str_interp("${dir}/match_score_matrix.csv"))){
            
            setwd(dir)
            
            df_msm = data.frame(fread("match_score_matrix.csv"))
            df_msm = df_msm %>% select(-V1)
            
            n_exceeding = sum(diag(as.matrix(df_msm)) > min_match_score)
            perc_exceeding = n_exceeding/nrow(df_msm)
            
            result = c(ratio, round(min_match_score, 3),
                       combination, 'all data',
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

results_df$`Panel size` = 192672

# TO BE FIXED
#clean_df = clean_up_table(results_df) 
#write.csv(clean_df, "ratios_and_percentages_TableS8.csv", row.names=FALSE)
