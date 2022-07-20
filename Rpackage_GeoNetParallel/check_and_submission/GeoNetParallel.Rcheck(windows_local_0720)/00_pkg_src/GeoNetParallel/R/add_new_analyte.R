#' Add new water chemistry data at new sampling locations
#' @description This function can be used to add new water chemistry data at new sampling locations and run the statistical inference for all the provided stream location. This function can only be used after “whole_pipe” is executed. Execution of this function will overwrite the previously generated output. Therefore, it is advised to preserve a copy of the old outputs in case you need it.
#' @param path_total This is the same root file path as that of the function “whole_pipe”. The previous output will be overwritten so backup is recommended.
#' @param df_analyte_to_append_filepath This is the full path to the csv file containing new water chemistry data for the given stream. The file formatting requirement is the same as that of “filename_chem” in “whole_pipe”.
#' @param upstream_thresh The same parameter as that in "whole_pipe".
#' @param downstream_lower_thresh The same parameter as that in "whole_pipe".
#' @param downstream_upper_thresh The same parameter as that in "whole_pipe".
#' @param num_cores The same parameter as that in "whole_pipe".
#' @param permission The same parameter as that in "whole_pipe".
#'
#' @return The same return files as those in "whole_pipe".
#' @note This function assumes “whole_pipe” function has already been run and the outputs generated from it lie in “path_total” folder. Please refer to our \href{https://github.com/HANDS-Research-Group/GeoNetParallel/tree/main/Rpackage_GeoNetParallel}{website} for examples.
#' @export
#'
#' @useDynLib GeoNetParallel
#' @importFrom Rcpp sourceCpp
#' @importFrom foreach %do%
#' @importFrom foreach %dopar%

add_new_analyte <- function(path_total, df_analyte_to_append_filepath,
                            upstream_thresh = 5,
                            downstream_lower_thresh = 0,
                            downstream_upper_thresh = 50,
                            num_cores = 1, permission = "no") {
  if (!(permission == "yes")){
    print("You need to explicit set the variable permission to make sure that you know what will happen. Please refer to the help page for function whole_pipe.")
    return (0)
  }
  #########handle the global variable issues###############
  df_anpoll_processed = 0
  df_analyte_nodeID_aggregated = 0
  igraph_river_whole = 0
  projected_nodeIDs_vec = 0
  polluter_vertex_ids = 0
  index = 0
  df_polluter_nodeID_aggregated = 0
  anpoll_edgelist = 0
  total_edgelist_character_modified = 0
  test_result = 0
  ########################################################
  file_path = path_total
  n_chunks = num_cores
  ##################################################
  df_analyte_to_append <- utils::read.csv(df_analyte_to_append_filepath)

  names(df_analyte_to_append) <- c("conc", "lon", "lat", "date")

  df_analyte_to_append$date <- as.Date(df_analyte_to_append$date, format = "%m/%d/%y")

  two_step_anpoll_mapper_list <- two_step_anpoll_mapper_update(df_anpoll_preprocessed = df_analyte_to_append,
                                                               output_path = file_path)

  load(file = paste0(file_path, "anpoll_files/df_anpoll_processed.RData"))
  load(file = paste0(file_path,"analyte_files/df_analyte_nodeID_aggregated.RData"))
  graph_subset_indicator <- F
  df_list_analyte_obj <- df_list_analyte_generator(df_analyte_processed = df_anpoll_processed[which(df_anpoll_processed$anpoll_indicator=="analyte"),],output_path = file_path)

  analyte_vertex_ids <- vertex_IDs_generator(df_anpoll_nodeID_aggregated = df_analyte_nodeID_aggregated,output_path = file_path,analyte=T,polluter=F,graph_subset=graph_subset_indicator)

  ########################################################################################################
  ## Creating shortest path edgelist by defining the segments and epochs and running in parallel 6 segments in each epoch with foreach::foreach

  ## Loading the analyte and polluter vertex ids
  load(file = paste0(file_path,"analyte_files/analyte_vertex_ids.RData"))
  load(file = paste0(file_path,"polluter_files/polluter_vertex_ids.RData"))
  ## Loading the vector of projected node IDs
  load(file = paste0(file_path,"polluter_files/projected_nodeIDs_vec.RData"))
  ## Loading the igraph object for whole river network
  load(file = paste0(file_path,"common_files_modified/igraph_river_whole.RData"))
  ## Getting the vertex IDs for the projected node IDs
  projected_vertex_ids <- which(igraph::V(igraph_river_whole)$name%in%projected_nodeIDs_vec)

  ## Combining the analyte, polluter and projected vertex ids
  analyte_polluter_projected_vertex_ids <- sort(unique(c(analyte_vertex_ids,polluter_vertex_ids,projected_vertex_ids)))

  ## Sourcing the modular functions for analysis

  ## Defining number of chunks to be parallelized, chunk size and indices inside each chunk
  chunk_size <- ceiling(length(analyte_polluter_projected_vertex_ids)/n_chunks)
  indices_all <- 1:length(analyte_polluter_projected_vertex_ids)
  indices_chunk_list <- split(x = indices_all, ceiling(seq_along(indices_all)/chunk_size))

  cl <- parallel::makeCluster(spec = n_chunks)
  doParallel::registerDoParallel(cl)
  shortest_path_anpoll_edgelist_chunk <- foreach::foreach (index = 1:n_chunks)%do% {

    shortest_path_edgelist_creator_parallelized_2(index=index, n_chunks=n_chunks,file_path = file_path)
  }
  parallel::stopCluster(cl)

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

  graph_subset_indicator <- F
  ## Loading the dataframe "df_analyte_nodeID_aggregated"
  load(file = paste0(file_path, "analyte_files/df_analyte_nodeID_aggregated.RData"))
  ## Loading the dataframe "df_polluter_nodeID_aggregated"
  load(file = paste0(file_path, "polluter_files/df_polluter_nodeID_aggregated.RData"))

  ## Loading the analyte-polluter network edgelist "anpoll_edgelist"
  load(file = paste0(file_path, "anpoll_files/anpoll_edgelist.RData"))
  ## Loading the analyte-polluter network edgelist "anpoll_edgelist"
  load(file = paste0(file_path, "anpoll_files/shortest_path_anpoll_edgelist.RData"))
  ## Loading the "total_edgelist_character_modified"
  load(file = paste0(file_path, "common_files_modified/total_edgelist_character_modified.RData"))

  ## Loading the dataframe "df_anpoll_processed"
  load(file = paste0(file_path, "anpoll_files/df_anpoll_processed.RData"))
  ## Loading the vector of projected node IDs
  load(file = paste0(file_path, "polluter_files/projected_nodeIDs_vec.RData"))

  ## Subsetting the df_polluter_processed from df_anpoll_processed
  df_polluter_processed <- df_anpoll_processed[which(df_anpoll_processed$anpoll_indicator == "polluter"),]

  ## Getting the flow dist to list
  flow_dist_to_list<-list()

  flow_dist_to_list  <- vector(mode="list", length=length( df_polluter_nodeID_aggregated))

  ## Actual Code
  ## cl <- parallel::makeCluster(2)
  ## doParallel::registerDoParallel(cl)
  flow_dist_to_list<-foreach::foreach (index = 1:nrow(df_polluter_nodeID_aggregated), .packages="Rcpp", .noexport = "getDistance")%do% {
    wrapper_flow_dist_cal(index, df_polluter= df_polluter_nodeID_aggregated,
                          anpoll_edgelist = anpoll_edgelist,
                          shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist,
                          total_edgelist_character_modified =
                            total_edgelist_character_modified,
                          from_indicator = F, to_indicator = T,
                          flow_type = "flow_to", file_path = file_path)
  }
  ## parallel::stopCluster(cl)

  names(flow_dist_to_list)<-df_polluter_nodeID_aggregated$nodeID
  print(" flow_dist_to_list is done")
  save(flow_dist_to_list,file = paste0(file_path,"inference/flow_dist_to_list.RData"))
  ############begin of part 3#################################
  load(file = paste0(file_path,"inference/flow_dist_from_list.RData"))
  ## Loading the "flow_dist_to_list"
  load(file = paste0(file_path,"inference/flow_dist_to_list.RData"))
  load(file = paste0(file_path,"inference/flow_dist_from_list_projected.RData"))

  flow_dist_polluter_projected_list <- flow_dist_polluter_projected_cal(file_path)

  load(file = paste0(file_path,"anpoll_files/anpoll_edgelist.RData"))
  load(file = paste0(file_path,"polluter_files/projected_nodeIDs_list.RData"))
  df_threshold_dist_km <- data.frame("polluter_intersection" = numeric(),"upstream" = numeric(), "downstream_lower" = numeric(),"downstream_upper" = numeric())

  df_threshold_dist_km[1,] <- c(50, upstream_thresh, downstream_lower_thresh, downstream_upper_thresh)

  for(j in 1:nrow(df_threshold_dist_km)) {
    ########################################################################################################

    load(file = paste0(file_path,"polluter_files/df_polluter_processed.RData"))
    dir.create(path = paste0(file_path,"inference/test_results"))
    if(n_chunks>nrow(df_polluter_processed)){
      n_chunks<-nrow(df_polluter_processed)
    }
    #chunk_size<-ceiling(12/n_chunks)
    chunk_size <- ceiling(nrow(df_polluter_processed)/n_chunks)
    #indices_all<-1:12
    indices_all <- 1:nrow(df_polluter_processed)
    indices_chunk_list <- split(x = indices_all, ceiling(seq_along(indices_all)/chunk_size))

    test_result_list <- list()

    ## index=1
    ## wrapper_polluter_test(n_chunks=n_chunks,file_path = file_path)
    cl <- parallel::makeCluster(spec = n_chunks)
    doParallel::registerDoParallel(cl)
    test_result_list <- foreach::foreach (index = 1:n_chunks)%dopar% {
      wrapper_polluter_test(index=index, n_chunks=n_chunks,file_path = file_path, j=j, df_threshold_dist_km = df_threshold_dist_km)
    }
    parallel::stopCluster(cl)

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
  #####begin of part 4########################
  df_threshold_dist_km <- data.frame("polluter_intersection" = numeric(),"upstream" = numeric(), "downstream_lower" = numeric(),"downstream_upper" = numeric())

  df_threshold_dist_km[1,] <- c(50, upstream_thresh, downstream_lower_thresh, downstream_upper_thresh)

  load(file = paste0(file_path,"polluter_files/df_polluter_processed.RData"))
  df_polluter_processed_appended  <- df_polluter_processed
  df_polluter_processed_appended$County  <-  NA
  ## df_polluter_processed_appended[nrow(df_polluter_processed_appended)+1,]  <- df_polluter_processed_appended[nrow(df_polluter_processed_appended),]
  for(j in 1:nrow(df_threshold_dist_km)) {
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
