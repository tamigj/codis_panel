#!/usr/bin/env Rscript

library(data.table)
library(vcfR)
library(stringr)
library(dplyr)
library(tidyr)


#--------------------#
# PROCESS ARGUMENTS  #
#--------------------#
args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
  stop("You must specify the STR", call.=FALSE)
} else if (length(args) != 0) {
  str = as.character(args[1])
}


#------------------#
# SOURCE VARIABLES #
#------------------#
# Source the shell script to set environment variables
system("source /scratch/groups/noahr/tami/codis_panel/config.sh")

# Retrieve the environment variables in R
dir_data_processed = Sys.getenv("DIR_DATA_PROCESSED")
dir_output_sumstats = Sys.getenv("DIR_OUTPUT_SUMSTATS")


#------------#
# FUNCTIONS  #
#------------#

# Calculate payseur D'avg and composite functions
calculate_payseur_d = function(snp, str, genotypes){

  haplotype_freqs_df = calculate_haplotype_frequencies(genotypes, str, snp)

  haplotypes_df = haplotype_freqs_df %>%

    # Calculate STR AF
    group_by(STR_allele) %>%
    mutate(STR_af = sum(Haplotype_freq)) %>%
    ungroup() %>%

    # Calculate SNP AF
    group_by(SNP_allele) %>%
    mutate(SNP_af = sum(Haplotype_freq)) %>%
    ungroup() %>%

    # Calculate D = Pab - Pa*Pb
    mutate(D = Haplotype_freq - SNP_af*STR_af) %>%

    rowwise() %>%

    # Calculate Dmax for negative and positive D
    mutate(Dmax_neg = min(SNP_af*STR_af, (1-SNP_af)*(1-STR_af)),
           Dmax_pos = min(SNP_af*(1-STR_af), STR_af*(1-SNP_af))) %>%

    # Calculate D'
    mutate(D_prime = ifelse(D < 0, D/Dmax_neg, D/Dmax_pos)) %>%

    # Calculate intermediate product
    mutate(product = STR_af*SNP_af*abs(D_prime))

  payseur_d = sum(haplotypes_df$product)

  return(payseur_d)

}

calculate_haplotype_frequencies = function(genotypes, str, snp){

  snp_and_str_df = genotypes[rownames(genotypes) %in% c(str, snp), ]

  haplotypes = c()

  for (i in 1:ncol(snp_and_str_df)){

    str_hap_alleles = unlist(str_split(snp_and_str_df[str, i], "\\|", 2))
    snp_hap_alleles = unlist(str_split(snp_and_str_df[snp, i], "\\|", 2))

    h1 = paste(str_hap_alleles[1], "-", snp_hap_alleles[1], sep='')
    h2 = paste(str_hap_alleles[2], "-", snp_hap_alleles[2], sep='')

    haplotypes = c(haplotypes, h1, h2)

  }

  # Calculate haplotype frequencies
  haplotype_freqs_df = data.frame(table(haplotypes))
  names(haplotype_freqs_df) = c('Haplotypes', 'Counts')

  haplotype_freqs_df$Haplotype_freq =
    haplotype_freqs_df$Counts/sum(haplotype_freqs_df$Counts)

  haplotype_freqs_df = haplotype_freqs_df %>%
    separate(Haplotypes, into=c('STR_allele', 'SNP_allele'), remove=T)

  # Add zero haplotype frequencies
  str_alleles = unique(haplotype_freqs_df$STR_allele)
  snp_alleles = unique(haplotype_freqs_df$SNP_allele)

  haplotype_freqs_df = add_missing_haplotypes(str_alleles, snp_alleles,
                                              haplotype_freqs_df)

  return(haplotype_freqs_df)

}

add_missing_haplotypes = function(str_alleles, snp_alleles, haplotype_table) {

  # Create a data frame to store the new rows
  missing_haplotypes <- expand.grid(STR_allele = str_alleles, SNP_allele = snp_alleles)

  # Convert alleles to character to ensure proper comparison
  missing_haplotypes$STR_allele <- as.character(missing_haplotypes$STR_allele)
  missing_haplotypes$SNP_allele <- as.character(missing_haplotypes$SNP_allele)

  # Merge with haplotype_table to find missing combinations
  merged <- merge(missing_haplotypes, haplotype_table, by = c("STR_allele", "SNP_allele"), all.x = TRUE)

  # Replace NA counts and frequencies with 0
  merged[is.na(merged$Counts), "Counts"] <- 0
  merged[is.na(merged$Haplotype_freq), "Haplotype_freq"] <- 0

  # Order the columns as in the original haplotype table
  merged <- merged[, names(haplotype_table)]

  return(merged)
}


#--------------------#
# LOAD & PROCESS VCF #
#--------------------#
setwd(dir_data_processed)

vcf_file = str_interp("${str}_withSTR_GT.vcf")
vcf = read.vcfR(vcf_file)

genotypes = extract.gt(vcf, element = "GT")


#--------------------------#
# CALCULARE PAYSEUR D'avg  #
#--------------------------#
snps = rownames(genotypes)[-which(rownames(genotypes) == str)]

snp_df = data.frame(matrix(nrow=0, ncol=2))
names(snp_df) = c('SNP', 'Payseur_d')

for (snp in snps){

  payseur_d = calculate_payseur_d(snp, str, genotypes)
  snp_df[nrow(snp_df)+1, ] = c(snp, payseur_d)

}


#------------#
# SAVE FILE  #
#------------#
setwd(dir_output_sumstats)
write.csv(snp_df, str_interp("${str}_payseur_ld.csv"), row.names=FALSE)
