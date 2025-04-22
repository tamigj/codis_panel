rm(list=ls())

library(dplyr)
library(ggplot2)
library(data.table)
library(stringr)

dir_out = "/scratch/groups/noahr/tami/codis_panel/output/payseur_ld"

codis_strs = c('CSF1PO', 'D10S1248', 'D12S391', 'D13S317',
               'D18S51', 'D19S433', 'D1S1656', 'D22S1045',
               'D2S1338', 'D2S441', 'D3S1358', 'D5S818',
               'D7S820', 'D8S1179', 'FGA', 'TH01', 'TPOX', 'vWA')

counter=0

for (d in c(0.3, 0.5, 0.7)){
  for (n_snps in c(25, 50, 75, 100)){
    for (i in 1:10){
      for (str in codis_strs){
        
        tmp_df = data.frame(fread(str_interp("${dir_out}/payseur_d_${d}_${n_snps}_${i}_${str}.ld")))
        
        tmp_df$payseur_d = d
        tmp_df$n_snps = n_snps
        tmp_df$snp_id = i
        tmp_df$str = str 
        
        tmp_df = tmp_df %>%
          select(payseur_d, n_snps, snp_id, str, R2)
          
        if (counter == 0){
          main_df = tmp_df
        } else {
          main_df = rbind(main_df, tmp_df)
        }
        
        counter=counter+1
      }
    }
  }
}

write.csv(main_df, str_interp("${dir_out}/payseur_ld_summary.csv"), row.names=FALSE)