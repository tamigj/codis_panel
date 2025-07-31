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

top_panels_df = data.frame(cbind(top_panels, top_panels_txt))
names(top_panels_df) = c('code', 'text')

ratios = c(10^0, 10^1, 10^2, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9, 10^10,
           10^11, 10^12, 10^13, 10^14, 10^15, 10^16, 10^17)


#------------#
# FUNCTIONS  #
#------------#
organize_results_in_table_all_data = function(results_df){
  
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
                     panel, 192672,
                     n_exceeding, round(perc_exceeding, 3))
          
          results_df[nrow(results_df)+1, ] = result
          
        }
      }
      
    }
  }
  
}

# setwd("/scratch/groups/noahr/tami/codis_panel/output/rm_summaries/0.75/")
# write.csv(results_df, "ratios_and_percentages_full_data.csv", row.names=FALSE)

clean_all_data_df = organize_results_in_table_all_data(results_df) 


#----------------------#
# MERGE AND SAVE DATA  #
#----------------------#
clean_df = rbind(clean_combinations_df, clean_all_data_df[-c(1:2), ])
write.csv(clean_df, "fraction_true_matches_TableS8.csv", row.names=FALSE)


#----------------#
# PLOT FOR NOAH  #
#----------------#
results_df$`Minimum match score` = 
  as.numeric(results_df$`Minimum match score`)

ggplot(results_df, aes(x=as.numeric(Percent))) +
  geom_histogram(bins=25) +
  facet_wrap(~`Minimum match score`, scales='free') +
  theme_light() + 
  xlab("Percent") + ylab("Count")


#---------------------#
# FIXED WIDTH x-AXIS  #
#---------------------#

# Calculate center point for each facet and set 0.1 width ranges
range_data <- results_df %>%
  group_by(`Minimum match score`) %>%
  summarise(
    center = (min(as.numeric(Percent), na.rm = TRUE) + max(as.numeric(Percent), na.rm = TRUE)) / 2,
    .groups = 'drop'
  ) %>%
  mutate(
    xlim_min = pmax(center - 0.05, 0),  # Don't go below 0
    xlim_max = pmin(center + 0.05, 1)   # Don't go above 1
  ) %>%
  mutate(
    # Adjust if we hit the boundaries
    xlim_min = ifelse(xlim_max == 1, 1 - 0.1, xlim_min),
    xlim_max = ifelse(xlim_min == 0, 0.1, xlim_max)
  )

# Create individual plots with exactly 0.1 width for each facet
plot_list <- lapply(unique(results_df$`Minimum match score`), function(score) {
  data_subset <- results_df[results_df$`Minimum match score` == score, ]
  limits <- range_data[range_data$`Minimum match score` == score, ]
  
  ggplot(data_subset, aes(x=as.numeric(Percent))) +
    geom_histogram(bins=length(unique(as.numeric(data_subset$Percent)))) +
    coord_cartesian(xlim = c(limits$xlim_min, limits$xlim_max)) +
    theme_light() + 
    xlab("Percent") + 
    ylab("Count") +
    ggtitle(paste("Min match score:", score)) +
    theme(
      plot.title = element_text(size = 10),
      axis.text.x = element_text(size = 8)
    ) +
    scale_x_continuous(breaks = c(limits$xlim_min, 
                                  (limits$xlim_min + limits$xlim_max)/2, 
                                  limits$xlim_max),
                       labels = function(x) sprintf("%.2f", x))
})

# Arrange in grid
library(gridExtra)
do.call(grid.arrange, c(plot_list, ncol = 5))

