rm(list=ls())

library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(cowplot)


#------------#
# VARIABLES. #
#------------#
cont_palette_3 = c('#E1ACAC', '#C75B7A', '#921A40')
combo_palette_3 = c('#92C7CF', '#088395', '#22577E')


#------------#
# FUNCTIONS  #
#------------#

# Data processing
make_summary_table = function(main_df, experiment_txt, fraction_num){
  
  df_summary = main_df %>%
    filter(experiment == experiment_txt,
           fraction == fraction_num) %>%
    pivot_longer(cols = c(one_to_one_median:needle_in_haystack_max), 
                 names_to = "temp",
                 values_to = "value") %>%
    separate(temp, into = c("Strategy", "Measure"), sep = "(?<=one_to_one|SNPquery|STRquery|needle_in_haystack)_") %>%
    mutate(Strategy = case_when(
      Strategy == "one_to_one" ~ "One-to-one",
      Strategy == "SNPquery" ~ "SNP query",
      Strategy == "STRquery" ~ "STR query",
      Strategy == "needle_in_haystack" ~ "Needle-in-haystack"
    )) %>%
    pivot_wider(names_from = Measure, values_from = value) %>%
    
    select(-n_replicates)
  
  df_summary$Strategy = factor(df_summary$Strategy,
                               levels = c("One-to-one", 'SNP query',
                                          'STR query', 'Needle-in-haystack'))
  
  # Drop filter column if it contains all NA's
  is_all_nas <- function(x) {
    all(is.na(x))
  }
  
  if (is_all_nas(df_summary$filter)) {
    df_summary <- df_summary %>% select(-filter)
  }
  
  # Factors of MAF/Pop-MAF
  if (experiment_txt == 'maf'){
    df_summary$filter = factor(df_summary$filter,
                               levels = c('1%', '5%', '10%'))
  } else if (experiment_txt == 'popmaf'){
    df_summary$filter = factor(df_summary$filter,
                               levels = c('0%', '1%', '5%', '10%'))
  }
  
  return(df_summary)
  
}

make_clean_summary_table = function(df_sum, pivot_longer=FALSE){
  
  df_sum_clean = df_sum %>%
    
    # Round to 3 decimal points 
    mutate(median = round(median, 3),
           min = round(min, 3),
           max = round(max, 3)) %>%
    
    # If it's 1, just have 1
    mutate(median = ifelse(median == 1, "1", median)) %>%
    
    # Construct ranges
    rowwise() %>%
    mutate(range = str_interp("[${min} - ${max}]")) %>%
    ungroup()
  
  # Keep only relevant columns
  if ("filter" %in% names(df_sum)){
    df_sum_clean = df_sum_clean %>%
      select(experiment, filter, n_total_snps, Strategy, median, range)
    
  } else {
    df_sum_clean = df_sum_clean %>%
      select(experiment, n_total_snps, Strategy, median, range) %>%
      mutate(filter = NA)
  }
  
  if (pivot_longer){
    
    df_pivot = df_sum_clean %>%
      pivot_wider(names_from = Strategy, values_from = c(median, range), 
                  names_sep = "_")
    
    # Rearrange and rename columns as desired
    df_sum_clean = df_pivot %>%
      
      mutate(filter = paste(experiment, filter, sep='-')) %>%
      
      select(filter, n_total_snps, 
             `median_One-to-one`, `range_One-to-one`,
             `median_SNP query`, `range_SNP query`,
             `median_STR query`, `range_STR query`,
             `median_Needle-in-haystack`, `range_Needle-in-haystack`) %>%
      
      rename( Filter = filter,
             `SNP panel size` = n_total_snps,
             `Median (One-to-one)` = `median_One-to-one`,
             `Range (One-to-one)` = `range_One-to-one`,
             `Median (SNP query)` = `median_SNP query`,
             `Range (SNP query)` = `range_SNP query`,
             `Median (STR query)` = `median_STR query`,
             `Range (STR query)` = `range_STR query`,
             `Median (Needle-in-haystack)` = `median_Needle-in-haystack`,
             `Range (Needle-in-haystack)` = `range_Needle-in-haystack`) %>%
      
      arrange(Filter)
    
  }
  
  return(df_sum_clean)
  
}

report_min_panel = function(df_sum, fraction){
  
  tmp = df_sum %>%
    filter(Strategy == 'Needle-in-haystack') %>%
    filter(median > 0.99) %>%
    filter(n_total_snps != 'all')
  
  if (nrow(tmp) != 0){
    tmp = tmp %>%
      filter(n_total_snps == min(n_total_snps))
    
    min_panel = tmp$n_total_snps
    print(str_interp("Minimum panel = ${min_panel} SNPs; ${fraction*100}% test."))
  } else {
    print(str_interp("No such panel."))
  }
  
  
}

process_filter_column <- function(df, col_name) {
  
  df$experiment[df$experiment == 'random_all_snps'] = 'random'
  
  # Individual filters ---------------------------------------
  df$filter[df$filter == 62500] = '0.0625Mb'
  df$filter[df$filter == 125000] = '0.125Mb'
  df$filter[df$filter == 250000] = '0.25Mb'

  df$filter[df$filter == 0] = '0%'
  df$filter[df$filter == 0.01] = '1%'
  df$filter[df$filter == 0.05] = '5%'
  df$filter[df$filter == 0.1] = '10%'
  
  # Combinations of filters -----------------------------------
  # Helper function to perform the replacement
  replace_values <- function(value) {
    if (value == "0") return("0%")
    if (value == "0.01") return("1%")
    if (value == "0.05") return("5%")
    if (value == "0.1") return("10%")
    if (value == "250000") return("0.25Mb")
    if (value == "125000") return("0.125Mb")
    if (value == "62500") return("0.0625Mb")
    return(value)
  }
  
  df[[col_name]] <- sapply(df[[col_name]], function(x) {
    if (is.na(x)) {
      return(NA)
    } else {
      # Determine the number of underscores
      num_underscores <- nchar(gsub("[^_]", "", x))
      
      if (num_underscores == 4) {
        # Split the string by underscores for 4 underscores
        split_values <- unlist(strsplit(x, "_"))
        names(split_values) <- c("MAF", "Pop-MAF", "Distance", "D\'avg", "R2")
        
      } else if (num_underscores == 3) {
        # Split the string by underscores for 3 underscores
        split_values <- unlist(strsplit(x, "_"))
        names(split_values) <- c("MAF", "Pop-MAF", "Distance", "R2")
        
      } else {
        return(x)
      }
      
      # Filter out NA values
      valid_values <- split_values[split_values != "NA"]
      
      # Concatenate names and values
      concatenated_values <- sapply(names(valid_values), function(name) {
        paste(name, replace_values(valid_values[name]))
      })
      
      # Return the concatenated string
      return(paste(concatenated_values, collapse = ", "))
    }
  }, USE.NAMES = FALSE)
  
  return(df)
}

select_combinations = function(df_sum, selection_criteria){
  
  if (selection_criteria == 'main_figure'){
    df_sum = df_sum %>%
      filter(!(grepl("D'avg ", filter)))
    
  } else if (selection_criteria == 'supplementary'){
    df_sum = df_sum %>%
      filter(grepl("D'avg ", filter))
  }
  
  return(df_sum)
}


# Plotting 
plot_rm_baseline = function(df_sum, title_txt){
  
  # Mark comparable combinations
  comparable_conditions = (df_sum %>%
                             filter(Strategy == 'Needle-in-haystack') %>%
                             filter(median > 0.99))$n_total_snps
  df_sum = df_sum %>%
    mutate(comparable = ifelse(n_total_snps %in% comparable_conditions,
                               "Comparable RM accuracies to all SNPs", 
                               "Worse RM accuracies than all SNPs"))
  
  ## Plot 
  p = ggplot(df_sum, aes(x=as.numeric(n_total_snps), y=median)) +
    geom_hline(yintercept=1, linetype='dashed', color='grey') +
    geom_pointrange(aes(ymin=min, ymax=max)) +
    # Add group=1 to ensure the line connects all points regardless of fill
    geom_line(aes(group=1)) +  
    
    geom_point(aes(fill=comparable), size=3, shape=21) +  
    facet_wrap(~Strategy, nrow=4) +
    theme_light() +
    theme(axis.text = element_text(size=16),
          axis.title = element_text(size=18),
          axis.text.x = element_text(angle=60, hjust=1),
          strip.text = element_text(size=18, color='White'),
          strip.background = element_rect(fill = "black"),
          plot.title = element_text(size=22, face='bold'),
          legend.position = 'none',
          legend.text = element_text(size=18),
          legend.title = element_blank()) + 
    
    scale_x_continuous(breaks = df_sum$n_total_snps) +
    scale_x_log10(labels = scales::comma_format(), 
                  breaks = df_sum$n_total_snps) +
    
    scale_fill_manual(values = c('black', 'white')) +
    
    ylim(c(0,1)) +
    xlab("SNP panel size") +
    ylab("\nRecord-matching accuracy\n")
  
  plot(p)
  
}

plot_rm = function(df_sum, title_txt, filter_txt){
  
  # Mark comparable combinations
  comparable_conditions = (df_sum %>%
                             filter(Strategy == 'Needle-in-haystack') %>%
                             filter(median > 0.99))$n_total_snps
  
  df_sum = df_sum %>%
    mutate(comparable = ifelse(n_total_snps %in% comparable_conditions,
                               "Comparable RM accuracies to all SNPs", 
                               "Worse RM accuracies than all SNPs"), 
           comparable = factor(comparable,
                               levels = c("Worse RM accuracies than all SNPs",
                                          "Comparable RM accuracies to all SNPs")))
  
  # Add the signs before filters
  if (title_txt %in% c('MAF', 'Pop-MAF', 'D_avg')){
    df_sum$filter = paste("≥", df_sum$filter, sep='')
    df_sum$filter = str_replace_all(df_sum$filter,
                                    "≥0%", ">0%")
    
  } else if (title_txt == 'Distance to CODIS STR'){
    df_sum$filter = paste("≤", df_sum$filter, sep='')
  }
  
  # Start with an empty plot and add the y=1 line first
  p = ggplot() +
    geom_hline(yintercept=1, linetype='dashed', color='grey')
  
  # Add color if there is filter!
  if ("filter" %in% names(df_sum)){
    
    if (grepl("Distance", title_txt)){
      df_sum = df_sum %>%
        mutate(filter = factor(filter, levels = c('≤0.25Mb', '≤0.125Mb', '≤0.0625Mb')))
      
    } else if (title_txt == "MAF"){
      df_sum = df_sum %>%
        mutate(filter = factor(filter, levels = c('≥1%', '≥5%', '≥10%')))
      
    } else {
      df_sum = df_sum %>%
        mutate(filter = factor(filter))
    }
    
    # Define colors for each filter
    n_filters <- length(unique(df_sum$filter))
    filter_colors <- cont_palette_3[1:n_filters]
    names(filter_colors) <- levels(df_sum$filter)
    
    # Split data by comparable status
    df_comparable <- df_sum %>% filter(comparable == "Comparable RM accuracies to all SNPs")
    df_worse <- df_sum %>% filter(comparable == "Worse RM accuracies than all SNPs")
    
    # Add lines using the complete dataset
    p = p + geom_line(data = df_sum, 
                      aes(x = as.numeric(n_total_snps), y = median, 
                          color = filter, group = filter))
    
    # Add error bars
    p = p + geom_pointrange(data = df_sum, 
                            aes(x = as.numeric(n_total_snps), y = median,
                                ymin = min, ymax = max,
                                color = filter))
    
    # Add filled points for "Comparable" data
    p = p + geom_point(data = df_comparable, 
                       aes(x = as.numeric(n_total_snps), y = median,
                           color = filter, shape = filter),
                       size = 3, fill = NA)  # NA fill means it inherits from color
    
    # Add white-filled points for "Worse" data
    p = p + geom_point(data = df_worse, 
                       aes(x = as.numeric(n_total_snps), y = median,
                           color = filter, shape = filter),
                       size = 3, fill = "white")
  } else {
    
    # For the case without filters
    p = p + geom_line(data = df_sum, 
                      aes(x = as.numeric(n_total_snps), y = median)) +
      geom_pointrange(data = df_sum,
                      aes(x = as.numeric(n_total_snps), y = median,
                          ymin = min, ymax = max)) +
      geom_point(data = df_sum,
                 aes(x = as.numeric(n_total_snps), y = median),
                 size = 3, fill = "white", shape = 21)
  }
  
  p = p +
    scale_shape_manual(values = c(21, 22, 23)) + # Filled circle, square, diamond
    
    facet_wrap(~Strategy, nrow=1) +
    theme_light() +
    theme(axis.text = element_text(size=16),
          axis.title = element_text(size=18),
          axis.text.x = element_text(angle=60, hjust=1),
          strip.text = element_text(size=18, color='White'),
          strip.background = element_rect(fill = "black"),
          plot.title = element_text(size=22, face='bold'),
          legend.position = 'bottom',
          legend.text = element_text(size=16),
          legend.title = element_text(size=16)) +
    scale_x_continuous(breaks = df_sum$n_total_snps) +
    ylim(c(0,1)) +
    xlab("SNP panel size") +
    ylab("\nRecord-matching\naccuracy\n") +
    guides(color = guide_legend(nrow = 1),
           shape = guide_legend(nrow = 1))
  
  if (title_txt == 'D_avg' & filter_txt == 'D_avg'){
    p = p + ggtitle(expression(D*"'"[avg])) +
      labs(color = expression(D*"'"[avg]), shape = expression(D*"'"[avg]))
  } else {
    p = p + 
      ggtitle(title_txt) +
      labs(color = filter_txt, shape = filter_txt)
  }
  
  # If there are filters, add color
  if ("filter" %in% names(df_sum)){
    p = p +
      scale_color_manual(values = filter_colors)
  } else {
    p = p +
      theme(legend.position = 'none')
  }
  
  plot(p)
  
}

plot_rm_combos = function(df_sum, filter_txt){
  
  df_sum = df_sum %>%
    mutate(Segment = NA) %>%
    mutate(Segment = ifelse(startsWith(filter, "MAF 1%"), 
                            'MAF≥1%', Segment),
           
           Segment = ifelse(startsWith(filter, "MAF 5%"), 
                            'MAF≥5%', Segment),
           
           Segment = ifelse(startsWith(filter, "MAF 10%"), 
                            'MAF≥10%', Segment),
           
           Segment = ifelse(startsWith(filter, "Pop-MAF 0%"),
                            'Pop-MAF>0%', Segment),
           
           Segment = ifelse(startsWith(filter, "Pop-MAF 1%"), 
                            'Pop-MAF≥1%', Segment),
           
           Segment = ifelse(startsWith(filter, "Pop-MAF 5%"), 
                            'Pop-MAF≥5%', Segment))
  
  df_sum$Segment = factor(df_sum$Segment,
                          levels = c('MAF≥1%',
                                     'MAF≥5%',
                                     'MAF≥10%',
                                     'Pop-MAF>0%',
                                     'Pop-MAF≥1%',
                                     'Pop-MAF≥5%'))
  
  # Mark comparable combinations
  comparable_conditions = (df_sum %>%
                             filter(Strategy == 'Needle-in-haystack') %>%
                             filter(median > 0.99))$n_total_snps
  
  df_sum = df_sum %>%
    mutate(comparable = ifelse(n_total_snps %in% comparable_conditions,
                               "Comparable RM accuracies to all SNPs", 
                               "Worse RM accuracies than all SNPs"), 
           comparable = factor(comparable,
                               levels = c("Worse RM accuracies than all SNPs",
                                          "Comparable RM accuracies to all SNPs")))
  
  # Add distance filter
  df_sum = df_sum %>%
    
    mutate(distance = NA) %>%
    
    mutate(distance = ifelse(grepl("0\\.0625Mb", filter),
                             '≤0.0625Mb', distance),
           
           distance = ifelse(grepl("0\\.125Mb", filter) & !grepl("0\\.0625", filter),
                             '≤0.125Mb', distance),
           
           distance = ifelse(grepl("0\\.25Mb", filter) & !grepl("0\\.125", filter) & !grepl("0\\.0625", filter),
                             '≤0.25Mb', distance)) %>%
    
    mutate(distance = factor(distance, 
                             levels = c('≤0.25Mb', '≤0.125Mb', '≤0.0625Mb'))) %>%
    
    arrange(factor(distance, levels = c('≤0.25Mb', '≤0.125Mb', '≤0.0625Mb')))
  
  # Start with an empty plot and add the y=1 line first
  p = ggplot() +
    geom_hline(yintercept=1, linetype='dashed', color='grey')
  
  # Define colors for each distance
  n_distances <- length(unique(df_sum$distance))
  distance_colors <- combo_palette_3[1:n_distances]
  names(distance_colors) <- levels(df_sum$distance)
  
  # Split data by comparable status
  df_comparable <- df_sum %>% filter(comparable == "Comparable RM accuracies to all SNPs")
  df_worse <- df_sum %>% filter(comparable == "Worse RM accuracies than all SNPs")
  
  # Add lines using the complete dataset
  p = p + geom_line(data = df_sum, 
                    aes(x = as.numeric(n_total_snps), y = median, 
                        color = distance, group = distance))
  
  # Add error bars
  p = p + geom_pointrange(data = df_sum, 
                          aes(x = as.numeric(n_total_snps), y = median,
                              ymin = min, ymax = max,
                              color = distance))
  
  # Add filled points for "Comparable" data
  p = p + geom_point(data = df_comparable, 
                     aes(x = as.numeric(n_total_snps), y = median,
                         color = distance, shape = distance),
                     size = 3, fill = NA)  # NA fill means it inherits from color
  
  # Add white-filled points for "Worse" data
  p = p + geom_point(data = df_worse, 
                     aes(x = as.numeric(n_total_snps), y = median,
                         color = distance, shape = distance),
                     size = 3, fill = "white")
  
  p = p +
    scale_shape_manual(values = c(21, 22, 23)) + # Filled circle, square, diamond
    
    facet_grid(Segment~Strategy) +
    theme_light() +
    theme(axis.text = element_text(size=18),
          axis.title = element_text(size=18, face='plain'),
          axis.text.x = element_text(angle=60, hjust=1),
          strip.text = element_text(size=18, color='White'),
          strip.background = element_rect(fill = "black"),
          plot.title = element_text(size=22, face='bold', hjust=0.5),
          legend.position = 'bottom',
          legend.text = element_text(size=18),
          legend.title = element_text(size=18),
          panel.spacing = unit(1, "lines")) +
    scale_x_continuous(breaks = unique(df_sum$n_total_snps)) + 
    scale_color_manual(values = distance_colors) +
    ylim(c(0,1)) +
    xlab("SNP panel size") +
    ylab("Record-matching accuracy\n") +
    labs(color = filter_txt, shape = filter_txt) +
    guides(color = guide_legend(nrow = 1),
           shape = guide_legend(nrow = 1))
  
  plot(p)
}

plot_rm_combos_supplement = function(df_sum){
  
  combinations = c("Distance 0.0625Mb, D'avg 0.3", "Distance 0.0625Mb, D'avg 0.5",
                   "Distance 0.0625Mb, D'avg 0.7",
                   "Distance 0.125Mb, D'avg 0.3", "Distance 0.125Mb, D'avg 0.5" ,
                   "Distance 0.125Mb, D'avg 0.7",
                   "Distance 0.25Mb, D'avg 0.3", "Distance 0.25Mb, D'avg 0.5", 
                   "Distance 0.25Mb, D'avg 0.7",
                   
                   "MAF 1%, D'avg 0.3", "MAF 1%, D'avg 0.5",
                   
                   "MAF 1%, Distance 0.0625Mb, D'avg 0.3",
                   "MAF 1%, Distance 0.125Mb, D'avg 0.3",
                   "MAF 1%, Distance 0.25Mb, D'avg 0.3",
                   
                   "MAF 1%, Distance 0.125Mb, D'avg 0.5", 
                   "MAF 1%, Distance 0.25Mb, D'avg 0.5")
  
  df_sum = df_sum %>%
    mutate(Segment = NA) %>%
    mutate(Segment = ifelse(filter %in% combinations[1:3],
                            'Dist≤0.0625Mb', Segment),
           
           Segment = ifelse(filter %in% combinations[4:6],
                            'Dist≤0.125Mb', Segment),
           
           Segment = ifelse(filter %in% combinations[7:9],
                            'Dist≤0.25Mb', Segment),
           
           Segment = ifelse(filter %in% combinations[10:11],
                            "MAF and D'avg", Segment),
           
           Segment = ifelse(filter %in% combinations[12:14],
                            "D'avg≥0.3", Segment),
           
           Segment = ifelse(filter %in% combinations[15:16],
                            "D'avg≥0.5", Segment))
  
  df_sum$Segment = factor(df_sum$Segment,
                          levels = c('Dist≤0.0625Mb',
                                     'Dist≤0.125Mb',
                                     'Dist≤0.25Mb',
                                     "MAF and D'avg",
                                     "D'avg≥0.3",
                                     "D'avg≥0.5"))
  
  # Mark comparable combinations
  comparable_conditions = (df_sum %>%
                             filter(Strategy == 'Needle-in-haystack') %>%
                             filter(median > 0.99))$n_total_snps
  
  df_sum = df_sum %>%
    mutate(comparable = ifelse(n_total_snps %in% comparable_conditions,
                               "Comparable RM accuracies to all SNPs", 
                               "Worse RM accuracies than all SNPs"), 
           comparable = factor(comparable,
                               levels = c("Worse RM accuracies than all SNPs",
                                          "Comparable RM accuracies to all SNPs")))
  
  # Add D'avg filter
  df_sum = df_sum %>%
    mutate(d = NA) %>%
    
    mutate(d = ifelse(grepl("D'avg 0.3", filter),
                      '≥0.3', d),
           
           d = ifelse(grepl("D'avg 0.5", filter),
                      '≥0.5', d),
           
           d = ifelse(grepl("D'avg 0.7", filter),
                      '≥0.7', d)) %>%
    
    mutate(d = factor(d))
  
  # Add distance filter 
  df_sum = df_sum %>%
    
    mutate(distance = NA) %>%
    
    mutate(distance = ifelse(grepl("0\\.0625Mb", filter),
                             '≤0.0625Mb', distance),
           
           distance = ifelse(grepl("0\\.125Mb", filter) & !grepl("0\\.0625", filter),
                             '≤0.125Mb', distance),
           
           distance = ifelse(grepl("0\\.25Mb", filter) & !grepl("0\\.125", filter) & !grepl("0\\.0625", filter),
                             '≤0.25Mb', distance)) %>%
    
    mutate(distance = factor(distance, 
                             levels = c('≤0.25Mb', '≤0.125Mb', '≤0.0625Mb')))
  
  # For splitting data by comparable status
  df_dist_comparable <- df_sum %>% 
    filter(grepl("Dist", Segment), comparable == "Comparable RM accuracies to all SNPs")
  df_dist_worse <- df_sum %>% 
    filter(grepl("Dist", Segment), comparable == "Worse RM accuracies than all SNPs")
  
  df_maf_comparable <- df_sum %>% 
    filter(Segment == "MAF and D'avg", comparable == "Comparable RM accuracies to all SNPs")
  df_maf_worse <- df_sum %>% 
    filter(Segment == "MAF and D'avg", comparable == "Worse RM accuracies than all SNPs")
  
  df_combo_comparable <- df_sum %>% 
    filter(grepl("avg≥", Segment), comparable == "Comparable RM accuracies to all SNPs")
  df_combo_worse <- df_sum %>% 
    filter(grepl("avg≥", Segment), comparable == "Worse RM accuracies than all SNPs")
  
  ## Plot (part 1) - Distance and D'avg
  pA = ggplot() +
    geom_hline(yintercept=1, linetype='dashed', color='grey') +
    
    # Add lines and error bars
    geom_line(data = df_sum %>% filter(grepl("Dist", Segment)),
              aes(x = as.numeric(n_total_snps), y = median, 
                  color = factor(d), group = factor(d))) +
    
    geom_pointrange(data = df_sum %>% filter(grepl("Dist", Segment)),
                    aes(x = as.numeric(n_total_snps), y = median,
                        ymin = min, ymax = max, color = factor(d))) +
    
    # Add filled points for "Comparable" data
    geom_point(data = df_dist_comparable, 
               aes(x = as.numeric(n_total_snps), y = median,
                   color = factor(d), shape = factor(d)),
               size = 3, fill = NA) +
    
    # Add white-filled points for "Worse" data
    geom_point(data = df_dist_worse, 
               aes(x = as.numeric(n_total_snps), y = median,
                   color = factor(d), shape = factor(d)),
               size = 3, fill = "white") +
    
    scale_shape_manual(values = c(21, 22, 23)) + # Filled circle, square, diamond
    
    facet_grid(Segment~Strategy) +
    theme_light() +
    theme(axis.text = element_text(size=16),
          axis.title = element_text(size=18, face='plain'),
          axis.text.x = element_text(angle=60, hjust=1),
          strip.text = element_text(size=18, color='White'),
          strip.background = element_rect(fill = "black"),
          plot.title = element_text(size=22, face='bold'),
          legend.position = 'bottom',
          legend.text = element_text(size=16),
          legend.title = element_text(size=16),
          panel.spacing = unit(1, "lines")) +
    scale_x_continuous(breaks = unique(df_sum$n_total_snps)) + 
    labs(x = "SNP panel size", color = expression(D*"'"[avg]), 
         shape = expression(D*"'"[avg])) +
    ylim(c(0,1)) +
    xlab("SNP panel size") +
    ylab("\nRecord-matching accuracy\n") +
    ggtitle(expression(bold(paste("Distance and ", D*"'"[avg])))) +
    guides(color = guide_legend(nrow = 1),
           shape = guide_legend(nrow = 1)) +
    scale_color_manual(values=combo_palette_3)
  
  ## Plot (part 2) - MAF and D'avg
  pB = ggplot() +
    geom_hline(yintercept=1, linetype='dashed', color='grey') +
    
    # Add lines and error bars
    geom_line(data = df_sum %>% filter(Segment == "MAF and D'avg"),
              aes(x = as.numeric(n_total_snps), y = median, 
                  color = factor(filter), group = factor(filter))) +
    
    geom_pointrange(data = df_sum %>% filter(Segment == "MAF and D'avg"),
                    aes(x = as.numeric(n_total_snps), y = median,
                        ymin = min, ymax = max, color = factor(filter))) +
    
    # Add filled points for "Comparable" data
    geom_point(data = df_maf_comparable, 
               aes(x = as.numeric(n_total_snps), y = median,
                   color = factor(filter), shape = factor(filter)),
               size = 3, fill = NA) +
    
    # Add white-filled points for "Worse" data
    geom_point(data = df_maf_worse, 
               aes(x = as.numeric(n_total_snps), y = median,
                   color = factor(filter), shape = factor(filter)),
               size = 3, fill = "white") +
    
    scale_shape_manual(values = c(21, 22)) + # Filled circle, square
    
    facet_wrap(~Strategy, nrow=1) +
    
    scale_color_brewer(
      palette='Paired',
      labels = c(
        "MAF 1%, D'avg 0.3" = expression(paste("MAF≥1%, ", D*"'"["avg"], "≥0.3")),
        "MAF 1%, D'avg 0.5" = expression(paste("MAF≥1%, ", D*"'"["avg"], "≥0.5"))
      )
    ) +
    
    theme_light() +
    theme(axis.text = element_text(size=16),
          axis.title = element_text(size=18, face='plain'),
          axis.text.x = element_text(angle=60, hjust=1),
          strip.text = element_text(size=18, color='White'),
          strip.background = element_rect(fill = "black"),
          plot.title = element_text(size=22, face='bold'),
          legend.position = 'bottom',
          legend.text = element_text(size=16),
          legend.title = element_blank(),
          panel.spacing = unit(1, "lines")) +
    scale_x_continuous(breaks = unique(df_sum$n_total_snps)) + 
    labs(x = "SNP panel size") +
    ylim(c(0,1)) +
    xlab("SNP panel size") +
    ylab("\nRecord-matching\naccuracy") +
    ggtitle(expression(bold(paste("MAF and ", D*"'"[avg])))) +
    
    guides(color = guide_legend(nrow = 1, override.aes = list(shape = c(21, 22))),
           shape = "none")
  
  ## Plot (part 3) - MAF≥1%, distance and D'avg
  pC = ggplot() +
    geom_hline(yintercept=1, linetype='dashed', color='grey') +
    
    # Add lines and error bars
    geom_line(data = df_sum %>% filter(grepl("avg≥", Segment)),
              aes(x = as.numeric(n_total_snps), y = median, 
                  color = factor(distance), group = factor(distance))) +
    
    geom_pointrange(data = df_sum %>% filter(grepl("avg≥", Segment)),
                    aes(x = as.numeric(n_total_snps), y = median,
                        ymin = min, ymax = max, color = factor(distance))) +
    
    # Add filled points for "Comparable" data
    geom_point(data = df_combo_comparable, 
               aes(x = as.numeric(n_total_snps), y = median,
                   color = factor(distance), shape = factor(distance)),
               size = 3, fill = NA) +
    
    # Add white-filled points for "Worse" data
    geom_point(data = df_combo_worse, 
               aes(x = as.numeric(n_total_snps), y = median,
                   color = factor(distance), shape = factor(distance)),
               size = 3, fill = "white") +
    
    scale_shape_manual(values = c(21, 22, 23)) + # Filled circle, square, diamond
    
    facet_grid(Segment~Strategy) +
    
    scale_color_manual(values = combo_palette_3) + 
    labs(color="Distance", shape="Distance") +
    
    theme_light() +
    theme(axis.text = element_text(size=16),
          axis.title = element_text(size=18, face='plain'),
          axis.text.x = element_text(angle=60, hjust=1),
          strip.text = element_text(size=18, color='White'),
          strip.background = element_rect(fill = "black"),
          plot.title = element_text(size=22, face='bold'),
          legend.position = 'bottom',
          legend.text = element_text(size=16),
          legend.title = element_text(size=16),
          panel.spacing = unit(1, "lines")) +
    scale_x_continuous(breaks = unique(df_sum$n_total_snps)) + 
    ylim(c(0,1)) +
    xlab("SNP panel size") +
    ylab("\nRecord-matching accuracy\n") +
    ggtitle(expression(bold(paste("MAF≥1%, distance and ", D*"'"[avg])))) +
    guides(color = guide_legend(nrow = 1, override.aes = list(shape = c(21, 22, 23))),
           shape = "none")
  
  return(plot_grid(pB, pA, pC, 
                   labels = c("A", "B", "C"), 
                   ncol = 1, 
                   rel_heights = c(0.22, 0.47, 0.33),
                   label_size = 24,          
                   label_fontface = "bold"))
}


#---------------#
# LOAD RM DATA  #
#---------------#
fraction = 0.75

setwd(str_interp('~/Desktop/codis_panel/rm_summaries/${fraction}'))
experiments = c('random', 'random_all_snps', 
                'maf', 'popmaf', 'distance', 'payseur_d',
                'combinations')
  
for (exp in experiments){
    
  path = str_interp("${exp}_rm_summary.csv")
  
  if (file.exists(path)){
    df = read.csv(path)
      
    if (exp %in% experiments[1:2]){
      df$filter = NA
    } else if (grepl('combination', exp)){
      df = df %>%
        rename('filter' = combination)
    }
      
    if (exp == experiments[1]){
      main_df = df 
    } else {
      main_df = rbind(main_df, df)
    }
  }
    
}
  
main_df = process_filter_column(main_df, 'filter')

  
#--------#
# PLOTS  #
#--------#
setwd('~/Desktop/codis_panel/output_plots/')
  
# 1. BASELINE  ------------------------------------------------------------------
df_sum = make_summary_table(main_df, 'random', fraction)
df_sum_clean = make_clean_summary_table(df_sum,
                                        pivot_longer = TRUE)
  
p = plot_rm_baseline(df_sum, 'Baseline: Random selection of SNPs')

png(str_interp('baseline_${fraction}.png'), width=700, height=800)
plot(p)
dev.off()
  
write.csv(df_sum_clean, 
          str_interp('baseline_${fraction}.csv'), 
          row.names=FALSE)
  
report_min_panel(df_sum, fraction)
  
  
# 2. VARIANT CHARACTERISTICS ----------------------------------------------------
tmp_df = main_df %>%
  filter(n_snps_per_str %in% c(25, 50, 75, 100))
  
df_sum = make_summary_table(tmp_df, 'random', fraction)
pA = plot_rm(df_sum, 'Random', 'NA')
report_min_panel(df_sum, fraction)
  
df_sum1 = make_summary_table(main_df, 'maf', fraction)
pB = plot_rm(df_sum1, 'MAF', 'MAF')
report_min_panel(df_sum1, fraction)
  
df_sum2 = make_summary_table(main_df, 'popmaf', fraction)
pC = plot_rm(df_sum2, 'Pop-MAF', 'Pop-MAF')
report_min_panel(df_sum2, fraction)
  
df_sum3 = make_summary_table(main_df, 'distance', fraction)
pD = plot_rm(df_sum3, 'Distance to CODIS STR', 'Distance to CODIS STR (in bp)')
report_min_panel(df_sum3, fraction)
  
df_sum4 = make_summary_table(main_df, 'payseur_d', fraction)
pE = plot_rm(df_sum4, 'D_avg', 'D_avg')
report_min_panel(df_sum4, fraction)
  
# Save plots
png(str_interp('variant_characteristics_all_${fraction}.png'), 
    width=950, height=1300)
plot_grid(NULL, pA, NULL, pB, pC, pD, pE, 
          labels = c("", "A", "", "B", "C", "D", "E"), 
          ncol = 1, 
          rel_heights = c(0.05, 0.85, 0.1, 1, 1, 1, 1),
          label_size = 24,          
          label_fontface = "bold",
          label_y=1.05) 
dev.off()
  
df_sum_clean = make_clean_summary_table(rbind(df_sum %>% mutate(filter=NA),
                                              df_sum1, df_sum2, 
                                              df_sum3, df_sum4), 
                                        pivot_longer = TRUE)
  
write.csv(df_sum_clean, 
          str_interp('variant_characteristics_${fraction}.csv'), 
          row.names=FALSE)
  
  
# 3. COMBINATIONS ---------------------------------------------------------------
df_sum_all = make_summary_table(main_df, 'combinations', fraction)
df_sum = select_combinations(df_sum_all, 'main_figure')

png(str_interp("combinations_all_${fraction}.png"), 
    width=1000, height=1000)
pA = plot_rm_combos(df_sum, 'Distance to CODIS STR (in bp)')
dev.off()
  
df_sum_clean = make_clean_summary_table(df_sum, pivot_longer = TRUE)
  
write.csv(df_sum_clean, 
          str_interp('combinations_main_${fraction}.csv'), 
          row.names=FALSE)
  
  
# 4. COMBINATIONS (supplementary) ----------------------------------------------
df_sum = select_combinations(df_sum_all, 'supplementary')
  
png(str_interp("combinations_supplement_${fraction}.png"), 
    width=1000, height=1300)
plot_rm_combos_supplement(df_sum)
dev.off()
  
df_sum_clean = make_clean_summary_table(df_sum, pivot_longer = TRUE)
  
write.csv(df_sum_clean, 
          str_interp('combinations_supplement_${fraction}.csv'), 
          row.names=FALSE)
  
