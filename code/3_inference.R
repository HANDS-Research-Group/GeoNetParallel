library(tidyverse)

## Update the path to the folder
file_path <- "C:/GeoNet/GeoNet_2021_packageDataset/"

## Sourcing the modular functions for analysis
source(file = paste0(file_path, "code/BaSu_network_v16_modular_functions.R"))

load(file = paste0(file_path,"inference/flow_dist_from_list.RData"))
## Loading the "flow_dist_to_list"
load(file = paste0(file_path,"inference/flow_dist_to_list.RData"))
load(file = paste0(file_path,"inference/flow_dist_from_list_projected.RData"))

flow_dist_polluter_projected_list <- flow_dist_polluter_projected_cal(file_path)

load(file = paste0(file_path,"anpoll_files/anpoll_edgelist.RData"))
load(file = paste0(file_path,"polluter_files/projected_nodeIDs_list.RData"))
## Defining the downstream_threshold_dist_km lower and upper
                                        #downstream_threshold_dist_km<-matrix(c(0,10,10,50),2,2,byrow=TRUE)
df_threshold_dist_km <- data.frame("polluter_intersection" = numeric(),"upstream" = numeric(), "downstream_lower" = numeric(),"downstream_upper" = numeric())

df_threshold_dist_km[1,] <- c(5, 5, 0, 10)
df_threshold_dist_km[2,] <- c(45, 5, 0, 50)

for(j in 1:2) {
  ########################################################################################################
  
  load(file = paste0(file_path,"polluter_files/df_polluter_processed.RData"))
  print(length(df_polluter_processed))
  dir.create(path = paste0(file_path,"inference/test_results"))
  n_chunks <- 2
  #chunk_size<-ceiling(12/n_chunks)
  chunk_size <- ceiling(nrow(df_polluter_processed)/n_chunks)
  #indices_all<-1:12
  indices_all <- 1:nrow(df_polluter_processed)
  indices_chunk_list <- split(x = indices_all, ceiling(seq_along(indices_all)/chunk_size))
  
  test_result_list <- list()
  
  ## index=1
  ## wrapper_polluter_test(n_chunks=n_chunks,file_path = file_path)
  cl <- makeCluster(spec = n_chunks)
  registerDoParallel(cl)
  test_result_list <- foreach (index = 1:n_chunks)%dopar% {
    wrapper_polluter_test(index=index, n_chunks=n_chunks,file_path = file_path)
  }
  stopCluster(cl)
  
  df_polluter_processed_appended  <- df_polluter_processed
  polluter_test_matrix <- data.frame(matrix(NA, nrow(df_polluter_processed), 15))
  upstream_downstream_obs_list<-list()
  for (i in 1:nrow(df_polluter_processed)){
      load(file = paste0(file_path,"inference/test_results/test_result_",i,".RData"))
      ## Getting mean, median and no. of observations for upstream and downstream in the first 6 columns
      polluter_test_matrix[i, c(1:6)] <- test_result[[3]][c(1:6)]
      ## Storing the p values and test indicator for upstream temporal test
      if(as.numeric(test_result[[1]][3]) == 1){
        polluter_test_matrix[i,7] <- as.numeric(test_result[[1]][1])
        polluter_test_matrix[i,8] <- as.numeric(test_result[[1]][2])
      }
      polluter_test_matrix[i,9] <- as.numeric(test_result[[1]][3])
      ## Storing the p values and test indicator for downstream temporal test
      if(as.numeric(test_result[[2]][3]) == 1){
        polluter_test_matrix[i,10] <- as.numeric(test_result[[2]][1])
        polluter_test_matrix[i,11] <- as.numeric(test_result[[2]][2])
      }
      polluter_test_matrix[i,12] <- as.numeric(test_result[[2]][3])
      ## Storing the p values and test indicator for final upstream vs. downstream spatio-temporal test
      if(as.numeric(test_result[[3]][9]) == 1){
        polluter_test_matrix[i,13] <- as.numeric(test_result[[3]][7])
        polluter_test_matrix[i,14] <- as.numeric(test_result[[3]][8])
      }
      polluter_test_matrix[i,15]<-as.numeric(test_result[[3]][9])
      
      # Lat/Lon information
      polluter_test_matrix[i,16]<-df_polluter_processed_appended$lon[i]
      polluter_test_matrix[i,17]<-df_polluter_processed_appended$lat[i]
      polluter_test_matrix[i,18]<-df_polluter_processed_appended$lon_mapped[i]
      polluter_test_matrix[i,19]<-df_polluter_processed_appended$lat_mapped[i]
      
      ## Initializing and storing the observations for upstream and downstream
      upstream_downstream_obs_list[[i]] <- list()
      upstream_downstream_obs_list[[i]][[1]] <- test_result[[4]]
      upstream_downstream_obs_list[[i]][[2]] <- test_result[[5]]
      #print(i)
      print(polluter_test_matrix[i,])
  
  }
  colnames(polluter_test_matrix) <- c("upstream mean","downstream mean","upstream median","downstream median","number of obs upstream","number of obs downstream", "t test (up) p value", "Wilcoxon test (up) p value","test result indicator (up)","t test (down) p value", "Wilcoxon test (down) p value","test result indicator (down)","t test (updown) p value", "Wilcoxon test (updown) p value","test result indicator (updown)", "lon", "lat", "lon_mapped", "lat_mapped")
  save(polluter_test_matrix, file = paste0(file_path, "inference/polluter_test_matrix_",df_threshold_dist_km[j,"downstream_upper"],".RData"))
  ## save(polluter_test_matrix, file = paste0(file_path, "inference/polluter_test_matrix_","50",".RData"))
}

