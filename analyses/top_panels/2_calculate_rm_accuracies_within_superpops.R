rm(list=ls())

#-----------------#
# LOAD LIBRARIES  #
#-----------------#
library(stringr)
library(dplyr)
library(data.table)
library(RecordMatching)


#------------#
# FUNCTIONS  #
#------------#
comp.match.acc.by.pop <- function(mat, AFR_ids, AMR_ids, EAS_ids, EUR_ids, SAS_ids) {
  
  if (dim(mat)[1] != dim(mat)[2]) {
    stop("The input match score matrix must be a square matrix.")
  }
  
  # Remove anyone with missing data at all loci
  mat <- mat[(rowSums(mat^2) != 0), (rowSums(mat^2) != 0)]
  
  # Create list of population ID lists
  pop_lists <- list(
    "AFR" = AFR_ids,
    "AMR" = AMR_ids, 
    "EAS" = EAS_ids,
    "EUR" = EUR_ids,
    "SAS" = SAS_ids
  )
  
  # Initialize results
  results <- data.frame()
  
  # Calculate accuracies for each population
  for (pop_name in names(pop_lists)) {
    pop_ids <- pop_lists[[pop_name]]
    
    # Get indices of these IDs in the matrix (only those that exist)
    pop_indices <- which(rownames(mat) %in% pop_ids)
    
    if (length(pop_indices) < 2) {
      next  # Skip if less than 2 individuals in population
    }
    
    # One-to-one matching (Hungarian algorithm on full matrix)
    matched.LSAP <- clue::solve_LSAP(mat - min(mat) + 1, maximum = TRUE)
    hungarian.acc <- mean(matched.LSAP[pop_indices] == pop_indices)
    
    # One-to-many matching - SNP query
    pickSTR.acc <- mean(apply(mat[, pop_indices], 2, which.max) == pop_indices)
    
    # One-to-many matching - STR query  
    pickSNP.acc <- mean(apply(mat[pop_indices, ], 1, which.max) == pop_indices)
    
    # Needle-in-haystack matching
    pop_mat <- mat[pop_indices, pop_indices]
    onematch.acc <- mean(true.greater.than.false(mat, pop_indices))
    
    # Add to results
    results <- rbind(results, data.frame(
      'Population' = pop_name,
      'n_individuals' = length(pop_indices),
      'one_to_one' = hungarian.acc,
      'SNPquery' = pickSTR.acc,
      'STRquery' = pickSNP.acc,
      'needle_in_haystack' = onematch.acc
    ))
  }
  
  return(results)
}

true.greater.than.false <- function(mat, indices = NULL){
  if (is.null(indices)) {
    indices <- 1:nrow(mat)
  }
  
  trues <- diag(mat)[indices]
  falses <- c(as.numeric(mat[lower.tri(mat, diag = FALSE)]),
              as.numeric(mat[upper.tri(mat, diag = FALSE)]))
  maxfalse <- max(falses)
  return(trues > maxfalse)
}

compute_rm_accuracies_in_pops = function(panel){
  
  count = 0
  
  for (i in 1:10){
    for (j in 1:10){
      
      dir = str_interp("/scratch/groups/noahr/tami/codis_panel/output/experiments/0.75/run_${panel}_${i}_${j}")
      
      # Process MSM if file exists 
      if (file.exists(dir) & file.exists(str_interp("${dir}/match_score_matrix.csv"))){
        
        setwd(dir)
        
        df_gp = data.frame(fread("CSF1PO.GP.FORMAT"))
        ids = names(df_gp)[-c(1,2)]
        
        df_msm = data.frame(fread("match_score_matrix.csv"))
        
        df_msm = df_msm %>% select(-V1)
        rownames(df_msm) = ids
        colnames(df_msm) = ids
        
        #results_all = RecordMatching::comp.match.acc(as.matrix(df_msm))
        results = comp.match.acc.by.pop(as.matrix(df_msm), 
                                        AFR_ids, AMR_ids, EAS_ids, EUR_ids, SAS_ids)
        
        if (count == 0){
          main_results = results
        } else {
          main_results = rbind(main_results, results)
        }
        
        count = count+1
      }
      
    }
  }
  
  print("Summarizing results across 100 replicates...")
  
  summarized_results = main_results %>%

    rowwise() %>%
    mutate(n_correct_one = n_individuals*one_to_one,
           n_correct_SNP = n_individuals*SNPquery,
           n_correct_STR = n_individuals*STRquery,
           n_correct_nih = n_individuals*needle_in_haystack) %>%
    ungroup() %>%
    
    group_by(Population) %>%
    summarize(mean_n_individuals = mean(n_individuals),
              
              mean_one_to_one = round(sum(n_correct_one)/sum(n_individuals), 3),
              min_one_to_one = round(min(one_to_one), 3),
              max_one_to_one = round(max(one_to_one), 3),
              
              mean_SNPquery = round(sum(n_correct_SNP)/sum(n_individuals), 3),
              min_SNPquery = round(min(SNPquery), 3),
              max_SNPquery = round(max(SNPquery), 3),
              
              mean_STRquery = round(sum(n_correct_STR)/sum(n_individuals), 3),
              min_STRquery = round(min(STRquery), 3),
              max_STRquery = round(max(STRquery), 3),
              
              mean_needle_in_haystack = round(sum(n_correct_nih)/sum(n_individuals), 3),
              min_needle_in_haystack = round(min(needle_in_haystack), 3),
              max_needle_in_haystack = round(max(needle_in_haystack), 3))

  return(summarized_results)
  
}

make_clean_table = function(main_results){
  
  summarized_results = main_results %>%
    mutate(range_one_to_one = paste("[", min_one_to_one, " - ", max_one_to_one, "]", sep=''),
           range_SNPquery = paste("[", min_SNPquery, " - ", max_SNPquery, "]", sep=''),
           range_STRquery = paste("[", min_STRquery, " - ", max_STRquery, "]", sep=''),
           range_needle_in_haystack = paste("[", min_needle_in_haystack, " - ", 
                                            max_needle_in_haystack, "]", sep='')) %>%
    
    mutate(panel_size = panel_size*18) %>%
    
    select(combination, panel_size,
           Population, mean_n_individuals,
           mean_one_to_one, range_one_to_one,
           mean_SNPquery, range_SNPquery,
           mean_STRquery, range_STRquery,
           mean_needle_in_haystack, range_needle_in_haystack)
  
  # Replace combination with text
  summarized_results$combination = 
    top_panels_df$text[match(summarized_results$combination, top_panels_df$code)]
  
  return(summarized_results)
}


#------------#
# VARIABLES  #
#------------#
dir_data_raw = '/scratch/groups/noahr/tami/codis_panel/data/raw/'
superpopulations = c('AFR', 'AMR', 'EUR', 'EAS', 'SAS')

top_panels = c('combinations_0.05_NA_125000_NA',
               'combinations_0.1_NA_125000_NA',
               'combinations_NA_0_62500_NA',
               'combinations_NA_0_125000_NA',
               'combinations_NA_0.01_62500_NA',
               'combinations_NA_0.05_125000_NA')

top_panels_txt = c("MAF≥5%, Distance≤0.125Mb",
                   "MAF≥10%, Distance≤0.125Mb",
                   "Pop-MAF≥0%, Distance≤0.0625Mb",
                   "Pop-MAF≥0%, Distance≤0.125Mb",
                   "Pop-MAF≥1%, Distance≤0.0625Mb",
                   "Pop-MAF≥5%, Distance≤0.125Mb")

top_panels_df = data.frame(cbind(top_panels, top_panels_txt))
names(top_panels_df) = c('code', 'text')

panel_sizes = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)


#----------------------------#
# LOAD SUPERPOPULATIONS IDs  #
#----------------------------#
for (pop in superpopulations){
  df_tmp = data.frame(fread((str_interp("${dir_data_raw}/${pop}_ids_plink.txt"))))
  ids = df_tmp$IID
  
  assign(str_interp("${pop}_ids"), ids)
}


#----------------------------------------------------#
# COMPUTE RM ACCURACIES IN SUPERPOPS FOR TOP PANELS  #
#----------------------------------------------------#
for (combination in top_panels){
  for (panel_size in panel_sizes){
    
    panel = str_interp("${combination}_${panel_size}")
    print(panel)
    
    results = compute_rm_accuracies_in_pops(panel)
  
    results$combination = combination
    results$panel_size = panel_size
    
    if (combination == top_panels[1] & panel_size == panel_sizes[1]){
      main_results = results
    } else {
      main_results = rbind(main_results, results)
    }
    
  }
}

clean_results_df = make_clean_table(main_results)


#------------#
# SAVE DATA  #
#------------#
setwd("/scratch/groups/noahr/tami/codis_panel/output/rm_summaries/0.75/")
write.csv(main_results, "top_panels_by_populations_rm.csv",
          row.names=FALSE)

write.csv(clean_results_df, "top_panels_by_populations_rm_summary_TableS7.csv")
