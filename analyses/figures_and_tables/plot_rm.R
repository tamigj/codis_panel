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
cont_palette_3 = c('#088395', '#D29F80', '#921A40')


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

process_filter_column = function(df, col_name) {
  
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
plot_rm_baseline_horizontal = function(df_sum, title_txt){
  
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
    geom_errorbar(aes(ymin=min, ymax=max), width=0.1) +
    # Add group=1 to ensure the line connects all points regardless of fill
    geom_line(aes(group=1)) +  
    
    geom_point(aes(fill=comparable), size=3, shape=21) +  
    facet_wrap(~Strategy, nrow=1) +
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
    
    #scale_x_continuous(breaks = c(180, 450, 1800, 9000,
    #                              18,000, 192672)) +
    scale_x_log10(labels = scales::comma_format(), 
                  breaks = c(180, 450, 1800, 9000, 192672)) +
    
    scale_fill_manual(values = c('black', 'white')) +
    
    ylim(c(0,1)) +
    xlab("SNP panel size") +
    ylab("Record-matching\naccuracy\n")
  
  plot(p)
  
}

plot_rm = function(df_sum, title_txt, filter_txt){
  
  # PRE-PROCESSING ----------------------------------------------
  
  # Fix filter names
  if (title_txt %in% c('MAF', 'Pop-MAF', 'D_avg')){
    df_sum$filter = paste("≥", df_sum$filter, sep='')
    df_sum$filter = str_replace_all(df_sum$filter,
                                    "≥0%", ">0%")
    
  } else if (title_txt == 'Distance to CODIS STR'){
    df_sum$filter = paste("≤", df_sum$filter, sep='')
  }
  
  # Fix filter order
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
  }
  
  if ("filter" %in% names(df_sum)){
    
    # COMPARABLE COMBINATIONS ------------------------------------
    df_sum = df_sum %>%
      mutate(candidate = ifelse(Strategy == "One-to-one" & median == 1 |
                                  Strategy == "SNP query" & median == 1 |
                                  Strategy == 'STR query' & median == 1 |
                                  Strategy == 'Needle-in-haystack' & median > 0.99,
                                1, 0)) %>%
      group_by(n_total_snps, filter) %>%
      mutate(n = sum(candidate))  %>%
      mutate(comparable = ifelse(n == 4,
                                 "Comparable RM accuracies to all SNPs", 
                                 "Worse RM accuracies than all SNPs"), 
             comparable = factor(comparable,
                                 levels = c("Worse RM accuracies than all SNPs",
                                            "Comparable RM accuracies to all SNPs")))
    
    
    df_comparable <- df_sum %>% 
      mutate(median = ifelse(comparable == "Comparable RM accuracies to all SNPs",
                             median, 5))
    
    df_worse <- df_sum %>% 
      filter(comparable == "Worse RM accuracies than all SNPs")
    
    # PLOT -------------------------------------------------------
    p = ggplot() +
      
      geom_hline(yintercept=1, linetype='dashed', color='grey') +
      facet_wrap(~Strategy, nrow=1) +
      
      # Add lines (all data)
      geom_line(data = df_sum, 
                aes(x = as.numeric(n_total_snps), y = median, 
                    color = filter, group = filter)) +
      
      # Add error bars (all data)
      geom_errorbar(data = df_sum, 
                    aes(x = as.numeric(n_total_snps), 
                        ymin = min, ymax = max,
                        color = filter), width=50) +
      
      geom_point(data = df_worse,
                 aes(x=n_total_snps, y=median,
                     shape=filter, color=filter), fill='white', size=4) +
      
      geom_point(data = df_comparable,
                 aes(x=n_total_snps, y=median,
                     shape=filter, color=filter, fill=filter), size=4) +
      
      scale_shape_manual(values = c(25, 22, 24)) + 
      scale_fill_manual(values = cont_palette_3) +
      scale_color_manual(values = cont_palette_3) +
      
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
      ylim(c(0,1.05)) +
      xlab("SNP panel size") +
      ylab("\nRecord-matching\naccuracy\n")
      
    # Add plot title and legend title
    if (title_txt == 'D_avg' & filter_txt == 'D_avg'){
      
      p = p + 
        ggtitle(expression(bold(D*"'"[avg]))) +
        guides(fill = "none", 
               shape = guide_legend(title = expression(D*"'"[avg])),
               color = guide_legend(title = expression(D*"'"[avg])))
      
    } else {
       p = p + 
         ggtitle(title_txt) +
         guides(fill = "none", 
                shape = guide_legend(title = filter_txt),
                color = guide_legend(title = filter_txt))
    }
  
  } else {
    
    # PLOT (BASELINE) --------------------------------------------
    p = ggplot() +
      
      geom_hline(yintercept=1, linetype='dashed', color='grey') +
      facet_wrap(~Strategy, nrow=1) +
      
      geom_line(data = df_sum, 
                      aes(x = as.numeric(n_total_snps), y = median)) +
        
      geom_errorbar(data = df_sum,
                    aes(x = as.numeric(n_total_snps), y = median,
                        ymin = min, ymax = max), width=50) +

      geom_point(data = df_sum,
                 aes(x = as.numeric(n_total_snps), y = median),
                 size = 4, fill = "white", shape = 21) +
      
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
      ylim(c(0,1.05)) +
      xlab("SNP panel size") +
      ylab("\nRecord-matching\naccuracy\n")
    
  }
  
  plot(p)
  
}

plot_rm_combos = function(df_sum, filter_txt){
  
  # PRE-PROCESSING ----------------------------------------------
  
  # Create Segment variable
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
  
  # COMPARABLE COMBINATIONS ------------------------------------
  df_sum = df_sum %>%
    mutate(candidate = ifelse(Strategy == "One-to-one" & median == 1 |
                                Strategy == "SNP query" & median == 1 |
                                Strategy == 'STR query' & median == 1 |
                                Strategy == 'Needle-in-haystack' & median > 0.99,
                              1, 0)) %>%
    group_by(n_total_snps, filter) %>%
    mutate(n = sum(candidate))  %>%
    mutate(comparable = ifelse(n == 4,
                               "Comparable RM accuracies to all SNPs", 
                               "Worse RM accuracies than all SNPs"), 
           comparable = factor(comparable,
                               levels = c("Worse RM accuracies than all SNPs",
                                          "Comparable RM accuracies to all SNPs")))
  
  # Create datasets for plotting
  df_comparable <- df_sum %>% 
    mutate(median = ifelse(comparable == "Comparable RM accuracies to all SNPs",
                           median, 5))
  
  df_worse <- df_sum %>% 
    filter(comparable == "Worse RM accuracies than all SNPs")
  
  # PLOT -------------------------------------------------------
  p = ggplot() +
    
    geom_hline(yintercept=1, linetype='dashed', color='grey') +
    facet_grid(Segment~Strategy) +
    
    # Add lines (all data)
    geom_line(data = df_sum, 
              aes(x = as.numeric(n_total_snps), y = median, 
                  color = distance, group = distance)) +
    
    # Add error bars (all data)
    geom_errorbar(data = df_sum, 
                  aes(x = as.numeric(n_total_snps), 
                      ymin = min, ymax = max,
                      color = distance), width=50) +
    
    # Add white-filled points for "Worse" data
    geom_point(data = df_worse,
               aes(x=n_total_snps, y=median,
                   shape=distance, color=distance, size=distance), fill='white') +
    
    # Add colored-filled points for "Comparable" data
    geom_point(data = df_comparable,
               aes(x=n_total_snps, y=median,
                   shape=distance, color=distance, fill=distance, size=distance)) +
    
    # Set scales - using cont_palette_3 to match plot_rm
    scale_shape_manual(values = c(21, 23, 22)) + 
    scale_size_manual(values = c(3, 5, 6)) +
    scale_fill_manual(values = cont_palette_3[1:3], guide = "none") +
    scale_color_manual(values = cont_palette_3[1:3]) +
    
    # Formatting
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
    ylim(c(0,1.05)) +
    xlab("SNP panel size") +
    ylab("Record-matching accuracy\n") +
    guides(fill = "none", 
           size = guide_legend(title = filter_txt),
           shape = guide_legend(title = filter_txt),
           color = guide_legend(title = filter_txt))
  
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
  
  # PRE-PROCESSING ----------------------------------------------
  
  # Create Segment variable
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
  
  # Add D'avg filter - Order for legend display
  df_sum = df_sum %>%
    mutate(d = NA) %>%
    
    mutate(d = ifelse(grepl("D'avg 0.3", filter),
                      '≥0.3', d),
           
           d = ifelse(grepl("D'avg 0.5", filter),
                      '≥0.5', d),
           
           d = ifelse(grepl("D'avg 0.7", filter),
                      '≥0.7', d)) %>%
    
    # Define factor levels for legend order
    mutate(d = factor(d, levels = c('≥0.3', '≥0.5', '≥0.7')))
  
  # Add distance filter
  df_sum = df_sum %>%
    
    mutate(distance = NA) %>%
    
    mutate(distance = ifelse(grepl("0\\.0625Mb", filter),
                             '≤0.0625Mb', distance),
           
           distance = ifelse(grepl("0\\.125Mb", filter) & !grepl("0\\.0625", filter),
                             '≤0.125Mb', distance),
           
           distance = ifelse(grepl("0\\.25Mb", filter) & !grepl("0\\.125", filter) & !grepl("0\\.0625", filter),
                             '≤0.25Mb', distance)) %>%
    
    # Define factor levels for the distance variable
    mutate(distance = factor(distance, 
                             levels = c('≤0.25Mb', '≤0.125Mb', '≤0.0625Mb')))
  
  # COMPARABLE COMBINATIONS ------------------------------------
  df_sum = df_sum %>%
    mutate(candidate = ifelse(Strategy == "One-to-one" & median == 1 |
                                Strategy == "SNP query" & median == 1 |
                                Strategy == 'STR query' & median == 1 |
                                Strategy == 'Needle-in-haystack' & median > 0.99,
                              1, 0)) %>%
    group_by(n_total_snps, filter) %>%
    mutate(n = sum(candidate))  %>%
    mutate(comparable = ifelse(n == 4,
                               "Comparable RM accuracies to all SNPs", 
                               "Worse RM accuracies than all SNPs"), 
           comparable = factor(comparable,
                               levels = c("Worse RM accuracies than all SNPs",
                                          "Comparable RM accuracies to all SNPs")))
  
  # Create datasets for plotting with proper sorting
  
  # Panel B - Distance and D'avg - Sort in reverse order for plotting
  df_dist <- df_sum %>% 
    filter(grepl("Dist", Segment)) %>%
    # Sort so that rows with ≥0.7 come first, then ≥0.5, then ≥0.3
    arrange(Strategy, Segment, n_total_snps, desc(as.numeric(gsub("≥", "", d))))
  
  df_dist_comparable <- df_dist %>% 
    mutate(median = ifelse(comparable == "Comparable RM accuracies to all SNPs",
                           median, 5))
  df_dist_worse <- df_dist %>% 
    filter(comparable == "Worse RM accuracies than all SNPs")
  
  # Panel A - MAF and D'avg
  df_maf <- df_sum %>% 
    filter(Segment == "MAF and D'avg") %>%
    # Make MAF 1%, D'avg 0.5 come before MAF 1%, D'avg 0.3 in dataframe
    arrange(Strategy, n_total_snps, desc(as.numeric(gsub("≥", "", gsub(".*D'avg ", "", filter)))))
  
  df_maf_comparable <- df_maf %>%
    mutate(median = ifelse(comparable == "Comparable RM accuracies to all SNPs",
                           median, 5))
  df_maf_worse <- df_maf %>%
    filter(comparable == "Worse RM accuracies than all SNPs")
  
  # Panel C - MAF≥1%, distance and D'avg
  df_combo <- df_sum %>% filter(grepl("avg≥", Segment))
  df_combo_comparable <- df_combo %>%
    mutate(median = ifelse(comparable == "Comparable RM accuracies to all SNPs",
                           median, 5))
  df_combo_worse <- df_combo %>%
    filter(comparable == "Worse RM accuracies than all SNPs")
  
  # Size mappings for d and distance - match the legend order
  d_size_map <- c('≥0.3' = 3, '≥0.5' = 4.5, '≥0.7' = 6)
  distance_size_map <- c('≤0.25Mb' = 3, '≤0.125Mb' = 4.5, '≤0.0625Mb' = 6)
  
  # For panel A, get unique filters in the order they appear in the data
  maf_filters <- unique(df_maf$filter)
  maf_size_map <- setNames(c(6, 4.5), maf_filters)
  
  # PLOTS ------------------------------------------------------
  
  # For panel A, use first and third colors from cont_palette_3
  maf_colors <- c(cont_palette_3[3], cont_palette_3[1])
  names(maf_colors) <- maf_filters
  
  ## Plot (part 1) - MAF and D'avg (Now A)
  pA = ggplot() +
    geom_hline(yintercept=1, linetype='dashed', color='grey') +
    facet_wrap(~Strategy, nrow=1) +
    
    # Add lines (all data)
    geom_line(data = df_maf, 
              aes(x = as.numeric(n_total_snps), y = median, 
                  color = filter, group = filter)) +
    
    # Add error bars (all data)
    geom_errorbar(data = df_maf, 
                  aes(x = as.numeric(n_total_snps), 
                      ymin = min, ymax = max,
                      color = filter), width=50) +
    
    # Add white-filled points for "Worse" data
    geom_point(data = df_maf_worse,
               aes(x=n_total_snps, y=median,
                   shape=filter, color=filter, size=filter), fill='white') +
    
    # Add colored-filled points for "Comparable" data
    geom_point(data = df_maf_comparable,
               aes(x=n_total_snps, y=median,
                   shape=filter, color=filter, fill=filter, size=filter)) +
    
    # Set scales
    scale_shape_manual(values = c(21, 23)) + 
    scale_size_manual(values = maf_size_map) +
    scale_color_manual(
      values = maf_colors,
      labels = c(
        "MAF 1%, D'avg 0.5" = expression(paste("MAF≥1%, ", D*"'"["avg"], "≥0.5")),
        "MAF 1%, D'avg 0.3" = expression(paste("MAF≥1%, ", D*"'"["avg"], "≥0.3"))
      )
    ) +
    scale_fill_manual(values = maf_colors, guide = "none") +
    
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
    scale_x_continuous(breaks = unique(df_maf$n_total_snps)) + 
    ylim(c(0,1.05)) +
    xlab("SNP panel size") +
    ylab("\nRecord-matching\naccuracy") +
    ggtitle(expression(bold(paste("MAF and ", D*"'"[avg])))) +
    guides(fill = "none", 
           shape = "none",
           size = "none",
           color = guide_legend(override.aes = list(shape = c(21, 23), size = c(6, 3))))
  
  ## Plot (part 2) - Distance and D'avg (Now B)
  pB = ggplot() +
    geom_hline(yintercept=1, linetype='dashed', color='grey') +
    facet_grid(Segment~Strategy) +
    
    # Add lines (all data)
    geom_line(data = df_dist, 
              aes(x = as.numeric(n_total_snps), y = median, 
                  color = d, group = d)) +
    
    # Add error bars (all data)
    geom_errorbar(data = df_dist, 
                  aes(x = as.numeric(n_total_snps), 
                      ymin = min, ymax = max,
                      color = d), width=50) +
    
    # Add white-filled points for "Worse" data
    geom_point(data = df_dist_worse,
               aes(x=n_total_snps, y=median,
                   shape=d, color=d, size=d), fill='white') +
    
    # Add colored-filled points for "Comparable" data
    geom_point(data = df_dist_comparable,
               aes(x=n_total_snps, y=median,
                   shape=d, color=d, fill=d, size=d)) +
    
    # Set scales with all three colors
    scale_shape_manual(values = c(21, 23, 22)) + 
    scale_size_manual(values = d_size_map) +
    scale_fill_manual(values = cont_palette_3, guide = "none") +
    scale_color_manual(values = cont_palette_3) +
    
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
    scale_x_continuous(breaks = unique(df_dist$n_total_snps)) + 
    ylim(c(0,1.05)) +
    xlab("SNP panel size") +
    ylab("\nRecord-matching accuracy\n") +
    ggtitle(expression(bold(paste("Distance and ", D*"'"[avg])))) +
    guides(fill = "none", 
           shape = "none",
           size = "none",
           color = guide_legend(title = expression(D*"'"[avg]), 
                                override.aes = list(shape = c(21, 23, 22), 
                                                    size = c(3, 4.5, 6))))
  
  ## Plot (part 3) - MAF≥1%, distance and D'avg
  pC = ggplot() +
    geom_hline(yintercept=1, linetype='dashed', color='grey') +
    facet_grid(Segment~Strategy) +
    
    # Add lines (all data)
    geom_line(data = df_combo, 
              aes(x = as.numeric(n_total_snps), y = median, 
                  color = distance, group = distance)) +
    
    # Add error bars (all data)
    geom_errorbar(data = df_combo, 
                  aes(x = as.numeric(n_total_snps), 
                      ymin = min, ymax = max,
                      color = distance), width=50) +
    
    # Add white-filled points for "Worse" data
    geom_point(data = df_combo_worse,
               aes(x=n_total_snps, y=median,
                   shape=distance, color=distance, size=distance), fill='white') +
    
    # Add colored-filled points for "Comparable" data
    geom_point(data = df_combo_comparable,
               aes(x=n_total_snps, y=median,
                   shape=distance, color=distance, fill=distance, size=distance)) +
    
    # Set scales
    scale_shape_manual(values = c(21, 23, 22)) + 
    scale_size_manual(values = distance_size_map) +
    scale_fill_manual(values = cont_palette_3, guide = "none") +
    scale_color_manual(values = cont_palette_3) +
    
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
    scale_x_continuous(breaks = unique(df_combo$n_total_snps)) + 
    ylim(c(0,1.05)) +
    xlab("SNP panel size") +
    ylab("\nRecord-matching accuracy\n") +
    ggtitle(expression(bold(paste("MAF≥1%, distance, and ", D*"'"[avg])))) +
    guides(fill = "none", 
           shape = "none",
           size = "none",
           color = guide_legend(title = "Distance", 
                                override.aes = list(shape = c(21, 23, 22), 
                                                    size = c(3, 4.5, 6))))
  
  return(plot_grid(pA, pB, pC, 
                   labels = c("A", "B", "C"), 
                   ncol = 1, 
                   rel_heights = c(0.22, 0.47, 0.33),
                   label_size = 24,          
                   label_fontface = "bold"))
}

plot_payseur_ld = function(df){
  
  df = df %>%
    filter(R2 > 0.1) %>%
    group_by(payseur_d, n_snps, snp_id) %>%
    mutate(mean_R2 = mean(R2)) %>%
    select(payseur_d, n_snps, snp_id, mean_R2) %>%
    mutate(n_snps = factor(n_snps, levels=c(25, 50, 75, 100)),
           payseur_d = factor(payseur_d)) %>%
    distinct()
  
  # Plot 
  p = ggplot(df,
             aes(x=factor(payseur_d), y=mean_R2)) +
    geom_boxplot(fill='#9CA986') +
    facet_wrap(~n_snps, nrow=1) +
    theme_light() +
    theme(axis.title = element_text(size=14),
          axis.text = element_text(size=12),
          plot.title = element_text(size=15, face='bold', hjust=0.5),
          strip.text = element_text(size=12, color='white'),
          strip.background = element_rect(fill = "black", color = "black")) +
    xlab(expression(paste("Minimum ", D*"'"[avg]))) +
    ylab(expression(paste("Mean r"^2, " between SNPs across 10 SNP sets")))
  
  return(p)
  
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
  
p = plot_rm_baseline_horizontal(df_sum, 'Baseline: Random selection of SNPs')

png(str_interp('baseline_${fraction}_horizontal.png'), width=1000, height=250)
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


# 5. LD among D'avg SNP sets (supplementary) -----------------------------------
df = read.csv(str_interp(
  '~/Desktop/codis_panel/rm_summaries/${fraction}/payseur_ld_summary.csv'))

png("Figure_S3.png", width=700, height=350)
plot_payseur_ld(df)
dev.off()
