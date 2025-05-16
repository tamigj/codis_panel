rm(list=ls())

library(dplyr)

#------------#
# FUNCTIONS  #
#------------#
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
# LOAD DATA  #
#------------#
setwd('~/Desktop/codis_panel/snp_info')
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


#------------------------------------------#
# HOW MANY SNPs REMAIN AFTER EACH FILTER?  #
#------------------------------------------#
setwd('~/Desktop/codis_panel/output_plots/')

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


summary_df$Characteristic = factor(summary_df$Characteristic,
                                   levels = c('MAF', 'Pop-MAF',
                                              'Distance to STR',
                                              'Payseur D'))

write.csv(summary_df, 'snps_remaining_after_filters_TableS5.csv', row.names=FALSE)
