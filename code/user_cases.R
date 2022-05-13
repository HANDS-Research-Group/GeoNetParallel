library(tidyverse)

## Update the path to the folder
file_path <- "C:/GeoNet/GeoNet_2021_packageDataset/"

## Sourcing the modular functions for analysis
source(file = paste0(file_path, "code/BaSu_network_v16_modular_functions.R"))

inference <- function() {
    source(file = paste0(file_path, "code/3_inference.R"))
    source(file = paste0(file_path, "code/4_fdr.R"))
}
# The following function takes in as input the inference output on all
# the polluter sites already present in GeoNet and a dataframe with
# new temporal polluter information and the threshold distance in km
# to produce new inference result
add_temporal_polluter <- function(df_polluter_to_append_file_path) {
    df_threshold_dist_km <- data.frame("polluter_intersection" = numeric(),"upstream" = numeric(), "downstream_lower" = numeric(),"downstream_upper" = numeric())

    df_threshold_dist_km[1,] <- c(5, 5, 0, 10)
    df_threshold_dist_km[2,] <- c(45, 5, 0, 50)

    df_polluter_raw<-read.csv(file = df_polluter_to_append_file_path)
    str(df_polluter_raw)

    ## Preprocessing the analyte dataframe
    df_polluter_preprocessed<-df_polluter_raw[,c(2,3)]
    df_polluter_preprocessed$date<-as.Date(df_polluter_raw$date,format="%m/%d/%Y")
        
    load(file=paste0(file_path,"polluter_files/df_polluter_processed.RData"))
    
    df_polluter_preprocessed$lon_mapped <- NA
    df_polluter_preprocessed$lat_mapped <- NA
    df_polluter_preprocessed$nodeID <- NA

    df_polluter_preprocessed$lon_mapped <- df_polluter_processed$lon_mapped[which(df_polluter_preprocessed$lon==df_polluter_processed$lon &
                                                                                    df_polluter_preprocessed$lat==df_polluter_processed$lat)]
    df_polluter_preprocessed$lat_mapped <- df_polluter_processed$lat_mapped[which(df_polluter_preprocessed$lat==df_polluter_processed$lat &
                                                                                    df_polluter_preprocessed$lon==df_polluter_processed$lon)]
    df_polluter_preprocessed$nodeID <- df_polluter_processed$nodeID[which(df_polluter_preprocessed$lon==df_polluter_processed$lon & 
                                                                            df_polluter_preprocessed$lat==df_polluter_processed$lat)]

    graph_subset_indicator <- F
    df_polluter_processed_appended  <- df_polluter_preprocessed
    for (j in 1:nrow(df_threshold_dist_km)){
        polluter_test_matrix_to_append <- data.frame(matrix(NA, nrow(df_polluter_preprocessed), 19))
        upstream_downstream_obs_list<-list()
        ## Loop to get test results for spills in rank_subgraph
        for(i in 1:nrow(df_polluter_preprocessed)){
            if(graph_subset_indicator){
                test_result <- polluter_test_dist_time(df_polluter = df_polluter_preprocessed, polluter_lon = df_polluter_preprocessed$lon_mapped[i], polluter_lat = df_polluter_preprocessed$lat_mapped[i], polluter_projected_dist_km = df_threshold_dist_km[j,"polluter_intersection"],upstream_threshold_dist_km = df_threshold_dist_km[j,"upstream"], downstream_threshold_lower_dist_km = df_threshold_dist_km[j,"downstream_lower"], downstream_threshold_upper_dist_km = df_threshold_dist_km[j,"downstream_upper"],date_start=as.Date("1920-01-01"), date_end=as.Date("2017-12-31"), spill_date = df_polluter_processed$date[i], file_path = file_path)
            }else{
                test_result <- polluter_test_dist_time(df_polluter = df_polluter_preprocessed, polluter_lon = df_polluter_preprocessed$lon_mapped[i], polluter_lat = df_polluter_preprocessed$lat_mapped[i], polluter_projected_dist_km = df_threshold_dist_km[j,"polluter_intersection"],upstream_threshold_dist_km = df_threshold_dist_km[j,"upstream"], downstream_threshold_lower_dist_km = df_threshold_dist_km[j,"downstream_lower"], downstream_threshold_upper_dist_km = df_threshold_dist_km[j,"downstream_upper"], date_start=as.Date("1920-01-01"), date_end=as.Date("2017-12-31"), spill_date = df_polluter_processed$date[i], file_path = file_path)
            }
            print(i)
                                        #assign("last.warning", NULL, envir = baseenv())
            ## Getting mean, median and no. of observations for upstream and downstream in the first 6 columns
            polluter_test_matrix_to_append[i, c(1:6)] <- test_result[[3]][c(1:6)]
            ## Storing the p values and test indicator for upstream temporal test
            if(as.numeric(test_result[[1]][3]) == 1){
                polluter_test_matrix_to_append[i,7] <- as.numeric(test_result[[1]][1])
                polluter_test_matrix_to_append[i,8] <- as.numeric(test_result[[1]][2])
            }
            polluter_test_matrix_to_append[i,9] <- as.numeric(test_result[[1]][3])
            ## Storing the p values and test indicator for downstream temporal test
            if(as.numeric(test_result[[2]][3]) == 1){
                polluter_test_matrix_to_append[i,10] <- as.numeric(test_result[[2]][1])
                polluter_test_matrix_to_append[i,11] <- as.numeric(test_result[[2]][2])
            }
            polluter_test_matrix_to_append[i,12] <- as.numeric(test_result[[2]][3])
            ## Storing the p values and test indicator for final upstream vs. downstream spatio-temporal test
            if(as.numeric(test_result[[3]][9]) == 1){
                polluter_test_matrix_to_append[i,13] <- as.numeric(test_result[[3]][7])
                polluter_test_matrix_to_append[i,14] <- as.numeric(test_result[[3]][8])
            }
            polluter_test_matrix_to_append[i,15]<-as.numeric(test_result[[3]][9])
            
                                        # Lat/Lon information
            polluter_test_matrix_to_append[i,16]<-df_polluter_processed_appended$lon[i]
            polluter_test_matrix_to_append[i,17]<-df_polluter_processed_appended$lat[i]
            polluter_test_matrix_to_append[i,18]<-df_polluter_processed_appended$lon_mapped[i]
            polluter_test_matrix_to_append[i,19]<-df_polluter_processed_appended$lat_mapped[i]
            
            ## Initializing and storing the observations for upstream and downstream
            upstream_downstream_obs_list[[i]] <- list()
            upstream_downstream_obs_list[[i]][[1]] <- test_result[[4]]
            upstream_downstream_obs_list[[i]][[2]] <- test_result[[5]]
                                        #print(i)
            print(polluter_test_matrix_to_append[i,])
        }
        
        colnames(polluter_test_matrix_to_append) <- c("upstream mean","downstream mean","upstream median","downstream median","number of obs upstream","number of obs downstream", "t test (up) p value", "Wilcoxon test (up) p value","test result indicator (up)","t test (down) p value", "Wilcoxon test (down) p value","test result indicator (down)","t test (updown) p value", "Wilcoxon test (updown) p value","test result indicator (updown)", "lon", "lat", "lon_mapped", "lat_mapped")
        load(file = paste0(file_path,"inference/polluter_test_matrix_",df_threshold_dist_km[j,"downstream_upper"],".RData"))
        polluter_test_matrix = rbind(polluter_test_matrix, polluter_test_matrix_to_append)
        save(polluter_test_matrix, file = paste0(file_path, "inference/polluter_test_matrix_",df_threshold_dist_km[j,"downstream_upper"],".RData"))
    }
    
    ####################
    ## FDR
    ####################
    
    df_polluter_processed_appended <- rbind(df_polluter_processed, df_polluter_processed_appended)
    df_polluter_processed_appended$County  <-  NA
    ## df_polluter_processed_appended[nrow(df_polluter_processed_appended)+1,]  <- df_polluter_processed_appended[nrow(df_polluter_processed_appended),]
    for(j in 1:2) {
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
    
}

add_temporal_analyte <- function(df_analyte_to_append_filepath) {
    load(file = paste0(file_path, "anpoll_files/df_anpoll_processed.RData"))
    df_analyte_to_append <- read.csv(df_analyte_to_append_filepath)

    names(df_analyte_to_append) <- c("conc", "lon", "lat", "date")
    df_analyte_to_append$date <- as.Date(df_analyte_to_append$date, format = "%m/%d/%y")

    df_analyte_to_append$anpoll_indicator <- "analyte"
    df_analyte_to_append$lon_mapped <- NA
    df_analyte_to_append$lat_mapped <- NA
    df_analyte_to_append$nodeID <- NA

    df_analyte_to_append$lon_mapped <- df_anpoll_processed$lon_mapped[which(df_anpoll_processed$lon==df_analyte_to_append$lon &
                                                                            df_anpoll_processed$lat==df_analyte_to_append$lat &
                                                                            df_anpoll_processed$date==df_analyte_to_append$date)[1]]
    df_analyte_to_append$lat_mapped <- df_anpoll_processed$lat_mapped[which(df_anpoll_processed$lon==df_analyte_to_append$lon &
                                                                            df_anpoll_processed$lat==df_analyte_to_append$lat &
                                                                            df_anpoll_processed$date==df_analyte_to_append$date)[1]]
    df_analyte_to_append$nodeID <- df_anpoll_processed$nodeID[which(df_anpoll_processed$lon==df_analyte_to_append$lon &
                                                                    df_anpoll_processed$lat==df_analyte_to_append$lat)[1]]

    df_anpoll_processed <- rbind(df_anpoll_processed, df_analyte_to_append)
    save(df_anpoll_processed, file = paste0(file_path, "anpoll_files/df_anpoll_processed.RData"))

    df_analyte_to_append <- df_anpoll_processed[which(df_anpoll_processed$anpoll_indicator=="analyte"),]
    df_list_analyte_obj <- df_list_analyte_generator(df_analyte_to_append, file_path)

    inference()
}

## Function takes the dataframe with new temporal analyte information
## and re-runs the inference for all polluter sites
add_new_analyte <- function(df_analyte_to_append_filepath) {
    load(file = paste0(file_path, "anpoll_files/df_anpoll_processed.RData"))
    df_analyte_to_append <- read.csv(df_analyte_to_append_filepath)

    names(df_analyte_to_append) <- c("conc", "lon", "lat", "date")

    df_analyte_to_append$date <- as.Date(df_analyte_to_append$date, format = "%m/%d/%y")

    two_step_anpoll_mapper_list <- two_step_anpoll_mapper_update(df_anpoll_preprocessed = df_analyte_to_append, 
                                                                output_path = file_path)
    load(file = paste0(file_path,"analyte_files/df_analyte_nodeID_aggregated.RData"))
    graph_subset_indicator <- F
    df_list_analyte_obj <- df_list_analyte_generator(df_analyte_processed = df_anpoll_processed[which(df_anpoll_processed$anpoll_indicator=="analyte"),],output_path = file_path)

    analyte_vertex_ids <- vertex_IDs_generator(df_anpoll_nodeID_aggregated = df_analyte_nodeID_aggregated,output_path = file_path,analyte=T,polluter=F,graph_subset=graph_subset_indicator)

    ########################################################################################################
    ## Creating shortest path edgelist by defining the segments and epochs and running in parallel 6 segments in each epoch with foreach

    ## Loading the analyte and polluter vertex ids
    load(file = paste0(file_path,"analyte_files/analyte_vertex_ids.RData"))
    load(file = paste0(file_path,"polluter_files/polluter_vertex_ids.RData"))
    ## Loading the vector of projected node IDs
    load(file = paste0(file_path,"polluter_files/projected_nodeIDs_vec.RData"))
    ## Loading the igraph object for whole river network
    load(file = paste0(file_path,"common_files_modified/igraph_river_whole.RData"))
    ## Getting the vertex IDs for the projected node IDs
    projected_vertex_ids <- which(V(igraph_river_whole)$name%in%projected_nodeIDs_vec)

    ## Combining the analyte, polluter and projected vertex ids
    analyte_polluter_projected_vertex_ids <- sort(unique(c(analyte_vertex_ids,polluter_vertex_ids,projected_vertex_ids)))
    
    ## Sourcing the modular functions for analysis
    
    ## Defining number of chunks to be parallelized, chunk size and indices inside each chunk
    n_chunks <- 2
    chunk_size <- ceiling(length(analyte_polluter_projected_vertex_ids)/n_chunks)
    indices_all <- 1:length(analyte_polluter_projected_vertex_ids)
    indices_chunk_list <- split(x = indices_all, ceiling(seq_along(indices_all)/chunk_size))

    #cl <- parallel::makeCluster(spec = n_chunks)
    #doParallel::registerDoParallel(cl)
    shortest_path_anpoll_edgelist_chunk <- foreach (index = 1:n_chunks)%do% {
      source(file = paste0(file_path, "code/BaSu_network_v16_modular_functions.R"))
      file_path <- "C:/GeoNet/GeoNet_2021_packageDataset/"
      shortest_path_edgelist_creator_parallelized_2(index=index, n_chunks=n_chunks,file_path = file_path)
    }
    #parallel::stopCluster(cl)

    shortest_path_anpoll_edgelist <- unlist(x = shortest_path_anpoll_edgelist_chunk,recursive = F)

    ## Removing duplicate paths
    shortest_path_anpoll_edgelist <- unique(shortest_path_anpoll_edgelist)

    ## Removing the length 1 paths
    shortest_path_length_vec <- sapply(X = shortest_path_anpoll_edgelist,FUN = length)
    if(length(which(shortest_path_length_vec == 1)) >= 1){
        shortest_path_anpoll_edgelist[which(shortest_path_length_vec == 1)] <- NULL
    }
    
    ## saving this path edgelist
    save(shortest_path_anpoll_edgelist, file = paste0(file_path,"anpoll_files/shortest_path_anpoll_edgelist.RData"))

    ## Creating the analyte polluter network
    anpoll_network_list <- anpoll_network_creator(output_path = file_path)

    head(anpoll_network_list[[1]])
    inference()
}

## The following function to add both analyte and polluter together at a time.
add_temporal_analyte_polluter <- function(df_analyte_to_append_filepath, df_polluter_to_append_file_path) {
    add_temporal_polluter(df_polluter_to_append)
    add_temporal_analyte(df_analyte_to_append)
}


df_polluter_to_append_file_path<- paste0(file_path,"data/Pulltion_Site_temporal.csv")
df_analyte_to_append_filepath <- paste0(file_path,"data/Water_chemistry_new.csv")

add_temporal_polluter(df_polluter_to_append_file_path)
