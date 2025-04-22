rm(list=ls())

library(dplyr)
library(ggplot2)
library(tidyr)

dir_out = "/scratch/groups/noahr/tami/codis_panel/output/payseur_ld"
setwd(dir_out)

#------------------------#
# LOAD AND PROCESS DATA  #
#------------------------#
df = read.csv('payseur_ld_summary.csv')

df = df %>%
  filter(R2 > 0.1) %>%
  group_by(payseur_d, n_snps, snp_id) %>%
  mutate(mean_R2 = mean(R2)) %>%
  
  select(payseur_d, n_snps, snp_id, mean_R2) %>%
  mutate(n_snps = factor(n_snps, levels=c(25, 50, 75, 100)),
         payseur_d = factor(payseur_d)) %>%
  
  distinct()


#--------#
# PLOTS  #
#--------#
ggplot(df,
       aes(x=factor(payseur_d), y=mean_R2)) +
  geom_boxplot(fill='#9CA986') +
  facet_wrap(~n_snps, nrow=1) +
  theme_light() +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        plot.title = element_text(size=15, face='bold', hjust=0.5),
        strip.text = element_text(size=12, color='white'),
        strip.background = element_rect(fill = "black", color = "black")) +
  xlab("Minimum Payseur D") +
  ylab("Mean R2 across SNPs\n(within an experiment)") +
  ggtitle("The higher the Payseur D threshold\n the higher the LD R2 between the selected SNPs")

