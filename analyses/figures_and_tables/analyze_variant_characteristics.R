## REMOVE LD-PRUNING

rm(list=ls())

library(dplyr)
library(ggplot2)
library(scales)

#install.packages("reshape2")
library(reshape2)
library(corrplot)
library(cowplot)
library(RColorBrewer)

setwd('~/Desktop/codis_panel/snp_info')


#------------#
# FUNCTIONS  #
#------------#
custom_label <- function(x) {
  ifelse(x %in% c(0, 1), as.character(x), sprintf("%.2f", x))
}

count_maf = function(df, t){
  
  tmp = df %>%
    filter(MAF >= t)
  
  return(nrow(tmp))
}

count_popmaf = function(df, t){
  
  if (t == 0){
    tmp = df %>%
      filter(AFR_MAF > t,
             AMR_MAF > t,
             EAS_MAF > t,
             EUR_MAF > t,
             SAS_MAF > t)
  } else {
    tmp = df %>%
      filter(AFR_MAF >= t,
             AMR_MAF >= t,
             EAS_MAF >= t,
             EUR_MAF >= t,
             SAS_MAF >= t)
  }
  
  return(nrow(tmp))
}

count_distance = function(df, t){
  
  tmp = df %>%
    filter(abs(dist_to_STR) <= t)
  
  return(nrow(tmp))
}

count_payseur_d = function(df, t){
  
  tmp = df %>%
    filter(Payseur_d >= t)
  
  return(nrow(tmp))
}


#------------#
# VARIABLES  #
#------------#
custom_colors = colorRampPalette(brewer.pal(11, "RdBu"))


#------------#
# LOAD DATA  #
#------------#
files = list.files(pattern = "*_snps.csv")

for (file in files){
 
  codis_str = sub("_snps.csv", "", file)
  
  tmp_df = read.csv(file)
  tmp_df$CODIS_STR = codis_str 
  
  if (file == files[1]){
    df = tmp_df
  } else {
    df = rbind(df, tmp_df)
  }
  
}


#------------#
# PAYSEUR D  #
#------------#
setwd('~/Desktop/codis_panel/output_snp_info')

p = ggplot(df, aes(x=Payseur_d)) +
  geom_histogram() +
  theme_light() +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        plot.title = element_text(size=22, face='bold')) +
  xlab("Payseur D") +
  ylab("Count")

png('payseur_d_histogram.png', width=450, height=350)
plot(p)
dev.off()


p = ggplot(df, aes(x=Payseur_d)) +
  geom_histogram() +
  facet_wrap(~CODIS_STR, nrow=3, ncol=6) +
  theme_light() +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=18, color='White'),
        strip.background = element_rect(fill = "black"),
        plot.title = element_text(size=22, face='bold')) +
  xlab("Payseur D") +
  ylab("Count") +
  scale_x_continuous(labels = custom_label)

png('payseur_d_histograms_by_str.png', width=1100, height=600)
plot(p)
dev.off()


#-------------------------#
# PAYSEUR D and distance  #
#-------------------------#

# Bin the data into intervals of 10,000 along the x-axis
df_binned <- df %>%
  mutate(bin = floor(dist_to_STR / 10000) * 10000) %>%  # Create bins of 10,000
  group_by(bin) %>%
  summarize(mean_payseur_d = mean(Payseur_d, na.rm = TRUE),
            sd_payseur_d = sd(Payseur_d, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(bin_center = bin + 5000)  # Calculate the center of each bin

# Plot the mean with error bars representing the standard deviation
p = ggplot(df_binned, aes(x = bin_center, y = mean_payseur_d)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_payseur_d - sd_payseur_d, ymax = mean_payseur_d + sd_payseur_d), width = 0) +
  theme_light() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 22, face = 'bold')) +
  xlab("Distance to STR (binned by 10,000bp)") +
  ylab("Mean Payseur D (with SD)")

png('payseur_d_vs_distance.png', width=600, height=450)
plot(p)
dev.off()


# Where are those with Payseur D of 1?
p = ggplot(df %>% filter(Payseur_d == 1), 
       aes(x = dist_to_STR)) +
  geom_histogram(bins=100) +
  theme_light() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 22, face = 'bold', hjust=0.5)) +
  xlab("Distance to STR") +
  ylab("Count") +
  ggtitle("Subsetted to SNPs with Payseur D = 1") +
  scale_x_continuous(labels = comma)

png('payseur_d_1_distance_distribution.png', width=450, height=350)
plot(p)
dev.off()


#----------------#
# CORRELOGRAMS   #
#----------------#
df_selected = df %>%
  select(MAF, 
         AFR_MAF, AMR_MAF, EAS_MAF, EUR_MAF, SAS_MAF, 
         Payseur_d, dist_to_STR) %>%
  mutate(dist_to_STR = abs(dist_to_STR))

# Compute Pearson correlation matrix
pearson_corr <- cor(df_selected, use = "complete.obs", method = "pearson")

# Compute Spearman correlation matrix
spearman_corr <- cor(df_selected, use = "complete.obs", method = "spearman")

# Pearson correlogram 
png('pearson_correlogram.png', width=800, height=550)
corrplot(pearson_corr, method = "color", type = "upper", 
         tl.cex = 1, title = "Pearson Correlation", 
         mar = c(0, 0, 2, 0), 
         addCoef.col = "black",  # Set all text to black
         number.cex = 1.2,  # Slightly larger text size for correlation values
         tl.col = "black",  # Set axis labels color to black
         diag = FALSE,  # Remove the diagonal
         col = custom_colors(200))  # Use the full RdBu color palette
dev.off()

# Plot Spearman correlogram with all black text
png('spearman_correlogram.png', width=800, height=550)
corrplot(spearman_corr, method = "color", type = "upper", 
         tl.cex = 1, title = "Spearman Correlation", 
         mar = c(0, 0, 2, 0), 
         addCoef.col = "black",  # Set all text to black
         number.cex = 1.2,  # Slightly larger text size for correlation values
         tl.col = "black",  # Set axis labels color to black
         diag = FALSE,  # Remove the diagonal
         col = custom_colors(200))  # Use the full RdBu color palette
dev.off()

#-------------------------------#
# MAKE SUPPLEMENTARY FIGURE 2   #
#-------------------------------#
pA = corrplot(pearson_corr, method = "color", type = "upper", 
              tl.cex = 1, title = "Pearson Correlation", 
              mar = c(0, 0, 2, 0), 
              addCoef.col = "black",  # Set all text to black
              number.cex = 1.2,  # Slightly larger text size for correlation values
              tl.col = "black",  # Set axis labels color to black
              diag = FALSE,  # Remove the diagonal
              col = custom_colors(200))  # Use the full RdBu color palette

pB = ggplot(df_binned, aes(x = bin_center, y = mean_payseur_d)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_payseur_d - sd_payseur_d, ymax = mean_payseur_d + sd_payseur_d), width = 0) +
  theme_light() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 22, face = 'bold')) +
  xlab("Distance to STR (binned by 10,000bp)") +
  ylab("Mean Payseur D (with SD)")

plot_grid(pA, pB,
          labels = c("A", "B"),
          ncol = 2, 
          label_size = 24,          
          label_fontface = "bold") 

#------------------------------------------#
# HOW MANY SNPs REMAIN AFTER EACH FILTER?  #
#------------------------------------------#
summary_df = data.frame(matrix(nrow=0, ncol=3))
names(summary_df) = c("Characteristic", "Filter", "N")

# MAF
summary_df[nrow(summary_df)+1, ] = c('MAF', 0.01, count_maf(df, 0.01))
summary_df[nrow(summary_df)+1, ] = c('MAF', 0.05, count_maf(df, 0.05))
summary_df[nrow(summary_df)+1, ] = c('MAF', 0.1, count_maf(df, 0.1))

# Pop-MAF
summary_df[nrow(summary_df)+1, ] = c('Pop-MAF', 0, 
                                     count_popmaf(df, 0))
summary_df[nrow(summary_df)+1, ] = c('Pop-MAF', 0.01, 
                                     count_popmaf(df, 0.01))
summary_df[nrow(summary_df)+1, ] = c('Pop-MAF', 0.05, 
                                     count_popmaf(df, 0.05))

# Distance to STR
summary_df[nrow(summary_df)+1, ] = c('Distance to STR',
                                     '0.5 Mb', count_distance(df, 250000))
summary_df[nrow(summary_df)+1, ] = c('Distance to STR', 
                                     '0.25 Mb', count_distance(df, 125000))
summary_df[nrow(summary_df)+1, ] = c('Distance to STR', 
                                     '0.125 Mb', count_distance(df, 62500))

# Payseur D
summary_df[nrow(summary_df)+1, ] = c('Payseur D', 0.3, 
                                     count_payseur_d(df, 0.3))
summary_df[nrow(summary_df)+1, ] = c('Payseur D', 0.5, 
                                     count_payseur_d(df, 0.5))
summary_df[nrow(summary_df)+1, ] = c('Payseur D', 0.7, 
                                     count_payseur_d(df, 0.7))

# LD pruning 
summary_df[nrow(summary_df)+1, ] = c('LD pruning', 0.3, 97736)
summary_df[nrow(summary_df)+1, ] = c('LD pruning', 0.5, 108655)
summary_df[nrow(summary_df)+1, ] = c('LD pruning', 0.7, 118302)

summary_df$Characteristic = factor(summary_df$Characteristic,
                                   levels = c('MAF', 'Pop-MAF',
                                              'Distance to STR',
                                              'Payseur D',
                                              'LD pruning'))

p = ggplot(summary_df, aes(x=Filter, y=as.numeric(N))) + 
  geom_bar(stat='identity') +
  facet_wrap(~Characteristic, scales='free_x', nrow=1) +
  theme_light() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(angle = 60, hjust=1),
        strip.text = element_text(size=18, color='White'),
        strip.background = element_rect(fill = "black"),
        plot.title = element_text(size = 22, face = 'bold', hjust=0.5)) +
  xlab("") +
  ylab("Count") +
  ggtitle("Number of SNPs remaining after filtering") +
  scale_y_continuous(labels = comma)

png('n_snps_remaining.png', width=900, height=450)
plot(p)
dev.off()

write.csv(summary_df, 'Table_S5.csv', row.names=FALSE)
