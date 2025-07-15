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
    mutate(n_correct = n_individuals*needle_in_haystack) %>%
    ungroup() %>%
    
    group_by(Population) %>%
    summarize(mean_n_individuals = mean(n_individuals),
              
              median_one_to_one = median(one_to_one),
              min_one_to_one = round(min(one_to_one), 3),
              max_one_to_one = round(max(one_to_one), 3),
              
              median_SNPquery = median(SNPquery),
              min_SNPquery = round(min(SNPquery), 3),
              max_SNPquery = round(max(SNPquery), 3),
              
              median_STRquery = median(STRquery),
              min_STRquery = round(min(STRquery), 3),
              max_STRquery = round(max(STRquery), 3),
              
              median_needle_in_haystack = median(needle_in_haystack),
              min_needle_in_haystack = round(min(needle_in_haystack), 3),
              max_needle_in_haystack = round(max(needle_in_haystack), 3),
              
              mean_needle_in_haystack = sum(n_correct)/sum(n_individuals))
  
  summarized_results = summarized_results %>%
    mutate(range_one_to_one = paste("[", min_one_to_one, " - ", max_one_to_one, "]", sep=''),
           range_SNPquery = paste("[", min_SNPquery, " - ", max_SNPquery, "]", sep=''),
           range_STRquery = paste("[", min_STRquery, " - ", max_STRquery, "]", sep=''),
           range_needle_in_haystack = paste("[", min_needle_in_haystack, " - ", 
                                            max_needle_in_haystack, "]", sep='')) %>%
    
    select(Population, mean_n_individuals,
           median_one_to_one, range_one_to_one,
           median_SNPquery, range_SNPquery,
           median_STRquery, range_STRquery,
           median_needle_in_haystack, range_needle_in_haystack,
           mean_needle_in_haystack)
  
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

panel_sizes = c(50, 60, 70, 80, 90, 100)


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


#---------------------#
# SAVE DATA AND PLOT  #
#---------------------#
setwd("/scratch/groups/noahr/tami/codis_panel/output/rm_summaries/0.75/")
write.csv(main_results, "top_panels_by_populations.csv",
          row.names=FALSE)

p = ggplot(main_results, 
       aes(x=panel_size, y=mean_needle_in_haystack, color=Population)) + 
  geom_hline(yintercept=0.99, color='black', linetype='dashed') +
  geom_point() + geom_smooth() + 
  facet_wrap(~combination) + theme_light()

png("top_panels_by_populations_mean_needle_in_haystack.png", width=800, height=400)
plot(p)
dev.off()

p = ggplot(main_results, 
           aes(x=panel_size, y=median_needle_in_haystack, color=Population)) + 
  geom_hline(yintercept=0.99, color='black', linetype='dashed') +
  geom_point() + geom_smooth() + 
  facet_wrap(~combination) + theme_light()

png("top_panels_by_populations_median_needle_in_haystack.png", width=800, height=400)
plot(p)
dev.off()
