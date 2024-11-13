rm(list=ls())

library(dplyr)
library(ggplot2)
library(data.table)
library(stringr)

dir_out = "/scratch/groups/noahr/tami/codis_panel/output/payseur_ld"

sum_df = data.frame(matrix(nrow=0, ncol=6))
names(sum_df) = c('payseur_d', 'n_snps', 'snp_n', 'mean_ld', 'sd_ld', 'n')

for (d in c(0.3, 0.5, 0.7)){
  for (n_snps in c(25, 50, 75, 100)){
    for (i in 1:10){
     
      tmp_df = data.frame(fread(str_interp("${dir_out}/payseur_d_${d}_${n_snps}_${i}_D7S820.ld")))
      
      mean_r2 = mean(tmp_df$R2)
      sd_r2 = sd(tmp_df$R2)
      n = nrow(tmp_df)
      
      sum_df[nrow(sum_df)+1, ] = c(d, n_snps, i, mean_r2, sd_r2, n)
      
    }
  }
}

write.csv(sum_df, str_interp("${dir_out}/payseur_ld_summary.csv"), row.names=FALSE)