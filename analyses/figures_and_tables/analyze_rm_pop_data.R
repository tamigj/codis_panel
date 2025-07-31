rm(list=ls())

library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(cowplot)

setwd("~/Desktop/codis_panel/figures_and_tables")


#------------#
# VARIABLES. #
#------------#
cont_palette_5 = c('#088395', '#D29F80', '#921A40', '#708A58', "#725CAD")

top_panels = c('combinations_0.05_NA_125000_NA',
               'combinations_0.1_NA_125000_NA',
               'combinations_NA_0_62500_NA',
               'combinations_NA_0_125000_NA',
               'combinations_NA_0.01_62500_NA',
               'combinations_NA_0.05_125000_NA')

top_panels_txt = c("MAF≥5%,\nDistance≤0.125Mb",
                   "MAF≥10%,\nDistance≤0.125Mb",
                   "Pop-MAF≥0%,\nDistance≤0.0625Mb",
                   "Pop-MAF≥0%,\nDistance≤0.125Mb",
                   "Pop-MAF≥1%,\nDistance≤0.0625Mb",
                   "Pop-MAF≥5%,\nDistance≤0.125Mb")

top_panels_df = data.frame(cbind(top_panels, top_panels_txt))
names(top_panels_df) = c('code', 'text')


#------------#
# FUNCTIONS  #
#------------#

# Data processing
make_summary_table = function(main_df){
  
  df_summary = main_df %>%
    select(-mean_n_individuals) %>%
    pivot_longer(cols = c(mean_one_to_one:max_needle_in_haystack), 
                 names_to = "temp",
                 values_to = "value") %>%
    separate(temp, into = c("Measure", "Strategy"), 
             sep = "_(?=one_to_one|SNPquery|STRquery|needle_in_haystack)") %>%
    mutate(Strategy = case_when(
      Strategy == "one_to_one" ~ "One-to-one",
      Strategy == "SNPquery" ~ "SNP query",
      Strategy == "STRquery" ~ "STR query",
      Strategy == "needle_in_haystack" ~ "Needle-in-haystack"
    )) %>%
    pivot_wider(names_from = Measure, values_from = value)
  
  df_summary$Strategy = factor(df_summary$Strategy,
                               levels = c("One-to-one", 'SNP query',
                                          'STR query', 'Needle-in-haystack'))
  
  return(df_summary)
  
}

# Ploting 
plot_rm_combos = function(df_sum, filter_txt){
  
  # PRE-PROCESSING ----------------------------------------------
  df_sum$Population = factor(df_sum$Population,
                             levels = c('AFR', 'AMR', 'EAS', 'EUR', 'SAS'))
  
  df_sum$combination =
    top_panels_df$text[match(df_sum$combination, top_panels_df$code)]
  
  df_sum = df_sum %>%
    rename(filter = 'combination') %>%
    mutate(n_total_snps = panel_size*18) %>%
    select(-panel_size)
  
  # COMPARABLE COMBINATIONS ------------------------------------
  df_sum = df_sum %>%
    mutate(candidate = ifelse(Strategy == "One-to-one" & mean == 1 |
                                Strategy == "SNP query" & mean == 1 |
                                Strategy == 'STR query' & mean == 1 |
                                Strategy == 'Needle-in-haystack' & mean > 0.99,
                              1, 0)) %>%
    group_by(Population, n_total_snps, filter) %>%
    mutate(n = sum(candidate))  %>%
    mutate(comparable = ifelse(n == 4,
                               "Comparable RM accuracies to all SNPs", 
                               "Worse RM accuracies than all SNPs"), 
           comparable = factor(comparable,
                               levels = c("Worse RM accuracies than all SNPs",
                                          "Comparable RM accuracies to all SNPs")))
  
  # Create datasets for plotting
  df_comparable <- df_sum %>% 
    mutate(mean = ifelse(comparable == "Comparable RM accuracies to all SNPs",
                         mean, 5))
  
  df_worse <- df_sum %>% 
    filter(comparable == "Worse RM accuracies than all SNPs")
  
  # PLOT -------------------------------------------------------
  p = ggplot() +
    
    geom_hline(yintercept=1, linetype='dashed', color='grey') +
    facet_grid(filter~Strategy) +
    
    # Add lines (all data)
    geom_line(data = df_sum, 
              aes(x = as.numeric(n_total_snps), y = mean, 
                  color = Population, group = Population)) +
    
    # Add error bars (all data)
    geom_errorbar(data = df_sum, 
                  aes(x = as.numeric(n_total_snps), 
                      ymin = min, ymax = max,
                      color = Population), width=50) +
    
    # Add white-filled points for "Worse" data
    geom_point(data = df_worse,
               aes(x=n_total_snps, y=mean,
                   shape=Population, color=Population, size=Population), fill='white') +
    
    # Add colored-filled points for "Comparable" data
    geom_point(data = df_comparable,
               aes(x=n_total_snps, y=mean,
                   shape=Population, color=Population, fill=Population, size=Population)) +
    
    # Set scales - using cont_palette_3 to match plot_rm
    scale_shape_manual(values = c(21, 22, 23, 24, 25)) + 
    scale_size_manual(values = c(4, 4.2, 4.5, 4.7, 5)) +
    scale_fill_manual(values = cont_palette_5, guide = "none") +
    scale_color_manual(values = cont_palette_5) +
    
    # Formatting
    theme_light() +
    theme(axis.text = element_text(size=18),
          axis.title = element_text(size=18, face='plain'),
          axis.text.x = element_text(angle=60, hjust=1),
          strip.text = element_text(size=15, color='White'),
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
           color = guide_legend(override.aes = list(shape = c(21, 22, 23, 24, 25), 
                                                    size = c(4, 4.2, 4.5, 4.7, 5))),
           shape = "none",
           size = "none") +
    labs(color = filter_txt)
  
  plot(p)
}


#---------------#
# LOAD RM DATA  #
#---------------#
fraction = 0.75

setwd(str_interp('~/Desktop/codis_panel/rm_summaries/${fraction}'))
df = read.csv('top_panels_by_populations_rm.csv')


#--------#
# PLOTS  #
#--------#
setwd('~/Desktop/codis_panel/output_plots/')

# PREPARE DATA -----
df_sum = make_summary_table(df)

# PLOT -------------
png(str_interp("top_panels_superpopulations_FigureS3.png"), 
    width=1000, height=1050)
plot_rm_combos(df_sum, 'Super-population')
dev.off()
