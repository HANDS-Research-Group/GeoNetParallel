
## Update the path to the folder
file_path <- "C:/GeoNet/GeoNet_2021_packageDataset/"

## Sourcing the modular functions for analysis
source(file = paste0(file_path, "code/BaSu_network_v16_modular_functions.R"))

df_threshold_dist_km <- data.frame("polluter_intersection" = numeric(), "upstream" = numeric(), "downstream_lower"=numeric(), "downstream_upper" = numeric())

df_threshold_dist_km[1,] <- c(5, 5, 0, 10)
df_threshold_dist_km[2,] <- c(45, 10, 0, 50)

load(file = paste0(file_path,"polluter_files/df_polluter_processed.RData"))
df_polluter_processed_appended  <- df_polluter_processed
df_polluter_processed_appended$County  <-  NA
## df_polluter_processed_appended[nrow(df_polluter_processed_appended)+1,]  <- df_polluter_processed_appended[nrow(df_polluter_processed_appended),]
for(j in 2:2) {
## for (j in 2:nrow(df_threshold_dist_km)){
    
  #####################################################
  ## Doing the fdr analysis for >=15 pvalues in each test
  load(file = paste0(file_path,"inference/polluter_test_matrix_",df_threshold_dist_km[j,"downstream_upper"],".RData"))
  df_polluter_test <- fdr_analysis_wrapper(polluter_test_matrix = polluter_test_matrix,alpha = 0.1,county=df_polluter_processed_appended$County, file_path = file_path)
  df_polluter_test_mean <- df_polluter_test[[1]]
  df_polluter_test_median <- df_polluter_test[[2]]
  save(df_polluter_test_mean, file = paste0(file_path,"inference/df_polluter_test_mean_", df_threshold_dist_km[j, "downstream_upper"], ".RData"))
  save(df_polluter_test_median, file = paste0(file_path,"inference/df_polluter_test_median_", df_threshold_dist_km[j, "downstream_upper"], ".RData"))
}

########################################################################################################
## Chunk parallezing the i for loop (as in above i.e. looping over different polluting events in df_polluter_processed_appended) for getting test results for different polluters using foreach

########################################################################################################
