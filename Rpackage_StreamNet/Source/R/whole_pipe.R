#' Some title
#' @description A description of this function with inline latex \eqn{2^3=8}
#' test for displayed equations \deqn{\frac{1}{2}\sum_{i=1}^n\log(x+y)}.
#' We can also use links, tables, lists and Character formattings such as bold and italic.
#' @param path_total description of this parameter
#' @param filename_chem description of this parameter
#' @param filename_poll description of this parameter
#' @param upstream_thresh description of this parameter
#' @param downstream_lower_thresh description of this parameter
#' @param downstream_upper_thresh description of this parameter
#' @param num_cores description of this parameter
#' @param permission description of this parameter
#'
#' @return description of the return value
#' @export
#'
#' @useDynLib Rpackage0621
#' @importFrom Rcpp sourceCpp
#' @importFrom foreach %do%
#' @importFrom foreach %dopar%

whole_pipe = function(path_total,
                      filename_chem = "Water_Chemistry_R_Package.csv",
                      filename_poll = "Pollution_Site_R_Package.csv",
                      upstream_thresh  = 5,
                      downstream_lower_thresh = 0,
                      downstream_upper_thresh = 50,
                      num_cores=1, permission = "no"){
  if (!(permission == "yes")){
    print("You need to explicit set this variable to make sure that you know what will happen.")
    return (0)
  }
  # upstream_thresh <- 5
  # downstream_lower_thresh <- 0
  # downstream_upper_thresh <- 50
  file_path = path_total
  n_chunks = num_cores
  ########################begin of part 1################################
  shape = rgdal::readOGR(dsn=paste0(file_path, "data/Flowline_R_Package"), layer="Flowline_R_Package")
  save(shape, file = paste0(file_path, "shape.RData"))
  ########################begin of part 2################################
  dir.create(path = paste0(file_path,"common_files"))
  dir.create(path = paste0(file_path,"common_files_modified"))
  dir.create(path = paste0(file_path,"analyte_files"))
  dir.create(path = paste0(file_path,"polluter_files"))
  dir.create(path = paste0(file_path,"anpoll_files"))
  dir.create(path = paste0(file_path,"anpoll_files/flow_from"))
  dir.create(path = paste0(file_path,"anpoll_files/flow_projected"))
  dir.create(path = paste0(file_path,"anpoll_files/flow_to"))
  dir.create(path = paste0(file_path,"anpoll_files/shortest_path_chunks"))
  dir.create(path = paste0(file_path,"anpoll_files/projected_missing_ids"))
  dir.create(path = paste0(file_path,"inference"))
  dir.create(path = paste0(file_path,"inference/test_results_j1"))
  dir.create(path = paste0(file_path,"inference/test_results_j2"))

  ###############################################################################################
  ## Loading the river net shape file
  load(file = paste0(file_path,"shape.RData"))

  ## Running the function "shapefile_pre_processing" to save the first three important files
  shapefile_pre_processing_list <- shapefile_pre_processing(shape_obj = shape,output_path = file_path)

  ###############################################################################################
  ## Loading the preprocessed analyte dataframe
  specCond <- utils::read.csv(file = paste0(file_path, 'data/', filename_chem))
  df_analyte_preprocessed = specCond[,c(1, 2, 3, 4)]
  names(df_analyte_preprocessed) <- c("conc", "lon", "lat", "date")

  df_analyte_preprocessed$date <- as.Date(df_analyte_preprocessed$date, format = "%m/%d/%y")


  ## analyte_data <- read.csv(file = paste0(file_path, 'data/Cl_data.csv'))
  ## df_analyte_preprocessed = analyte_data[,c(1, 2, 3)]
  ## df_analyte_preprocessed$date = as.Date(analyte_data$Date, "%m/%d/%y")
  ## names(df_analyte_preprocessed) <- c("conc", "lon", "lat", "date")
  ## save(df_analyte_preprocessed, file=paste0(file_path, "analyte_files/df_analyte_preprocessed.RData"))

  ###################################################
  ## Preprocessing the polluter dataframe
  df_polluter_raw<-utils::read.csv(file = paste0(file_path,"data/",filename_poll))
  names(df_polluter_raw) <- c("county", "lon", "lat", "date")
  ## Preprocessing the analyte dataframe
  df_polluter_preprocessed<-df_polluter_raw[,c(2,3)]
  df_polluter_preprocessed$date<-as.Date(df_polluter_raw$date,format="%m/%d/%Y")
  ## saving this dataframe
  save(df_polluter_preprocessed,file = paste0(file_path,"polluter_files/df_polluter_preprocessed.RData"))

  ####################################################
  ## Combining the analyte and polluter dataframes
  df_anpoll_preprocessed <- data.frame(matrix(NA,nrow(df_analyte_preprocessed)+nrow(df_polluter_preprocessed),4))
  df_anpoll_preprocessed <- data.frame(df_analyte_preprocessed,"anpoll_indicator"=rep("analyte",nrow(df_analyte_preprocessed)),stringsAsFactors = F)
  df_anpoll_preprocessed[(nrow(df_analyte_preprocessed)+1):(nrow(df_analyte_preprocessed)+nrow(df_polluter_preprocessed)),c(2,3,4)] <- df_polluter_preprocessed
  df_anpoll_preprocessed[(nrow(df_analyte_preprocessed)+1):(nrow(df_analyte_preprocessed)+nrow(df_polluter_preprocessed)),"anpoll_indicator"] <- "polluter"


  ## tallying the number of polluter sources
  # df_anpoll_preprocessed %>% as_tibble() %>% group_by(anpoll_indicator) %>% summarise(n = n())
  #  summarise(group_by(as_tibble(df_anpoll_preprocessed), anpoll_indicator), n())
  ## saving this dataframe
  save(df_anpoll_preprocessed, file = paste0(file_path,"anpoll_files/df_anpoll_preprocessed.RData"))

  ########################################################################################################
  ## Loading the dataframe "df_anpoll_preprocessed"
  load(file = paste0(file_path,"anpoll_files/df_anpoll_preprocessed.RData"))

  ## Performing the two step dynamic mapping procedure that can map all C-PP locations to the base river map
  two_step_anpoll_mapper_list <- two_step_anpoll_mapper(df_anpoll_preprocessed = df_anpoll_preprocessed,output_path = file_path)

  ## Calculating the path distances of all streams inside "stream_list_modified"
  stream_path_dist_vec <- stream_path_dist_cal(output_path = file_path)

  ## Loading the dataframe "df_anpoll_processed"
  load(file = paste0(file_path,"anpoll_files/df_anpoll_processed.RData"))

  ## Getting the 2nd MOST IMPORTANT dataframe and other lists for analyte
  df_list_analyte_obj <- df_list_analyte_generator(df_analyte_processed = df_anpoll_processed[which(df_anpoll_processed$anpoll_indicator=="analyte"),],output_path = file_path)

  ## Getting the 3rd MOST IMPORTANT dataframe for polluter
  df_polluter_nodeID_aggregated <- df_polluter_generator(df_polluter_processed = df_anpoll_processed[which(df_anpoll_processed$anpoll_indicator=="polluter"), c("nodeID","lon_mapped","lat_mapped","date")], output_path = file_path)

  ########################################################################################################
  ## Subsetting the river net and getting 1) "igraph_river_whole" 2) "graph order" and 3) "igraph_river_decomposed_list"
  river_net_subset_list <- river_net_subset(output_path = file_path)

  ## Defining graph_subset_indicator and rank_subgraph
  graph_subset_indicator <- F
  #rank_subgraph<-2

  ## Loading the df_analyte_nodeID_aggregated and df_polluter_nodeID_aggregated
  load(file = paste0(file_path,"analyte_files/df_analyte_nodeID_aggregated.RData"))
  load(file = paste0(file_path,"polluter_files/df_polluter_nodeID_aggregated.RData"))

  ########################################################################################################
  ########################################################################################################
  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###
  # {## Test code!!
  #   ## checking the output files of two step anpoll mapper function and others
  #   load(file=paste0(file_path,"anpoll_files/dist_mapped.RData"))
  #   str(dist_mapped)
  #   summary(dist_mapped)
  #   hist(dist_mapped)
  #   hist(dist_mapped[dist_mapped<200])
  #
  #   load(file=paste0(file_path,"anpoll_files/df_anpoll_processed.RData"))
  #   str(df_anpoll_processed)
  #   head(df_anpoll_processed)
  #   length(unique(df_anpoll_processed$nodeID))
  #   geosphere::distm(x = as.matrix(df_anpoll_processed[1,c("lon","lat")]),y = as.matrix(df_anpoll_processed[1,c("lon_mapped","lat_mapped")]))
  #
  #   load(file=paste0(file_path,"polluter_files/df_polluter_processed.RData"))
  #   str(df_polluter_processed)
  #   summary(tail(dist_mapped))
  #
  #   load(file=paste0(file_path,"common_files_modified/df_node_latlong_modified.RData"))
  #   str(df_node_latlong_modified)
  #   head(df_node_latlong_modified)
  #
  #   load(file=paste0(file_path,"common_files/df_node_latlong.RData"))
  #   str(df_node_latlong)
  #   length(unique(df_node_latlong$nodeID))
  #
  #   load(file = paste0(file_path,"common_files/stream_list.RData"))
  #   str(stream_list)
  #
  #   load(file = paste0(file_path,"common_files_modified/stream_list_modified.RData"))
  #   str(stream_list_modified)
  #
  #   load(file = paste0(file_path,"common_files_modified/stream_path_dist_vec.RData"))
  #   str(stream_path_dist_vec)
  #   summary(stream_path_dist_vec)
  #
  #   load(file = paste0(file_path,"analyte_files/df_analyte_nodeID_aggregated.RData"))
  #   str(df_analyte_nodeID_aggregated)
  #
  #   load(file = paste0(file_path,"analyte_files/list_analyte_time.RData"))
  #   str(list_analyte_time, max.level = 1)
  #
  #   load(file = paste0(file_path,"analyte_files/listID_nodeID_matrix.RData"))
  #   str(listID_nodeID_matrix)
  #   head(listID_nodeID_matrix)
  #
  #   load(file = paste0(file_path,"polluter_files/df_polluter_nodeID_aggregated.RData"))
  #   str(df_polluter_nodeID_aggregated)
  #
  #   load(file = paste0(file_path,"common_files_modified/igraph_river_decomposed_list.RData"))
  #   str(igraph_river_decomposed_list, max.level = 1)}

  ##  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###
  #######################################################################################################

  ## Getting the analyte_vertex_ids
  analyte_vertex_ids <- vertex_IDs_generator(df_anpoll_nodeID_aggregated = df_analyte_nodeID_aggregated,output_path = file_path,analyte=T,polluter=F,graph_subset=graph_subset_indicator)

  ## Getting the polluter_vertex_ids
  polluter_vertex_ids <- vertex_IDs_generator(df_anpoll_nodeID_aggregated = df_polluter_nodeID_aggregated,output_path = file_path,analyte=F,polluter=T, graph_subset=graph_subset_indicator)

  ## Getting the projected_nodeIDs_list
  projected_nodeIDs_list <- projected_nodeIDs_list_generator(file_path=file_path, projected_threshold_dist_km = 50)

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

  ## Defining number of chunks to be parallelized, chunk size and indices inside each chunk
  chunk_size <- ceiling(length(analyte_polluter_projected_vertex_ids)/n_chunks)
  indices_all <- 1:length(analyte_polluter_projected_vertex_ids)
  indices_chunk_list <- split(x = indices_all, ceiling(seq_along(indices_all)/chunk_size))

  cl <- parallel::makeCluster(spec = n_chunks)
  doParallel::registerDoParallel(cl)
  shortest_path_anpoll_edgelist_chunk <- foreach::foreach (index = 1:n_chunks)%dopar% {
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

  ########################################################################################################
  ####################################################
  ## Creating the analyte polluter network
  anpoll_network_list <- anpoll_network_creator(output_path = file_path)


  ######################################################################################################
  ######################################################################################################
  ############################################# Inference ##############################################
  ######################################################################################################
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
  nrow(df_polluter_processed)

  ## Subsetting the dataframe "df_polluter_nodeID_aggregated" for the
  ## case when river network is subsetted geographically, calculating
  ## the flow distance for each polluter from and to all connecting
  ## nodes in the anpoll network and storing them in a list
  {if(graph_subset_indicator){
    df_polluter_nodeID_aggregated_rank_subgaph<-df_polluter_nodeID_aggregated[which(df_polluter_nodeID_aggregated$rank_subgraph),]
    row.names(df_polluter_nodeID_aggregated_rank_subgaph)<-NULL
    ## Getting the flow dist from list
    flow_dist_from_list<-list()
    cl <- parallel::makeCluster(n_chunks)
    doParallel::registerDoParallel(cl)
    flow_dist_from_list<-foreach::foreach (index = 1:nrow(df_polluter_nodeID_aggregated_rank_subgaph))%dopar% {
      wrapper_flow_dist_cal(index=index, df_polluter=df_polluter_nodeID_aggregated_rank_subgaph,anpoll_edgelist = anpoll_edgelist, shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist, total_edgelist_character_modified = total_edgelist_character_modified, from_indicator = T, to_indicator = F, file_path = file_path)
    }
    parallel::stopCluster(cl)
    names(flow_dist_from_list)<-df_polluter_nodeID_aggregated_rank_subgaph$nodeID
    save(flow_dist_from_list,file = paste0(file_path,"inference/flow_dist_from_list.RData"))
    print("from_list is done")
    ## Getting the flow dist to list
    flow_dist_to_list<-list()
    cl <- parallel::makeCluster(n_chunks)
    doParallel::registerDoParallel(cl)
    flow_dist_to_list<-foreach::foreach (index = 1:nrow(df_polluter_nodeID_aggregated_rank_subgaph))%dopar% {
      wrapper_flow_dist_cal(index=index, df_polluter=df_polluter_nodeID_aggregated_rank_subgaph,anpoll_edgelist = anpoll_edgelist, shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist, total_edgelist_character_modified = total_edgelist_character_modified, from_indicator = F, to_indicator = T, file_path = file_path)
    }
    parallel::stopCluster(cl)
    names(flow_dist_to_list)<-df_polluter_nodeID_aggregated_rank_subgaph$nodeID
    save(flow_dist_to_list,file = paste0(file_path,"inference/flow_dist_to_list.RData"))
  }else{
    ## Getting the flow dist from list for polluters
    flow_dist_from_list_polluters<-list()

    ## ptm <- proc.time()
    print("From_list_polluter started")

    ## Uncomment the following section
    ## cl <- parallel::makeCluster(2)
    ## doParallel::registerDoParallel(cl)
    ##flow_dist_from_list_polluters <- foreach::foreach (index = 1:nrow(df_polluter_nodeID_aggregated)) %dopar%{
    flow_dist_from_list_polluters <- foreach::foreach (index = 1:nrow(df_polluter_nodeID_aggregated), .packages="Rcpp", .noexport = "getDistance") %do% {
      wrapper_flow_dist_cal(index, df_polluter = df_polluter_nodeID_aggregated,
                            anpoll_edgelist = anpoll_edgelist,
                            shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist,
                            total_edgelist_character_modified =
                              total_edgelist_character_modified,
                            from_indicator = T, to_indicator = F, flow_type = "flow_from",
                            file_path = file_path)
    }
    ## parallel::stopCluster(cl)
    names(flow_dist_from_list_polluters)<-df_polluter_nodeID_aggregated$nodeID
    save(flow_dist_from_list_polluters,
         file = paste0(file_path,"inference/flow_dist_from_list_polluters.RData"))
    print("from_list_polluter is done")


    projected_nodeIDs_vec = unique(projected_nodeIDs_vec)
    ## Getting the flow dist from list for projected node IDs of polluters
    df_projected_nodeIDs<-data.frame("nodeID"=projected_nodeIDs_vec,stringsAsFactors = F)
    flow_dist_from_list_projected<-list()

    ## cl <- parallel::makeCluster(2)
    ## doParallel::registerDoParallel(cl)
    flow_dist_from_list_projected<-foreach::foreach (index =1:nrow(df_projected_nodeIDs), .packages="Rcpp", .noexport = "getDistance")%do% {
      wrapper_flow_dist_cal(index, df_polluter=df_projected_nodeIDs,
                            anpoll_edgelist = anpoll_edgelist,
                            shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist,
                            total_edgelist_character_modified = total_edgelist_character_modified,
                            from_indicator = T, to_indicator = F,
                            flow_type = "flow_projected",
                            file_path = file_path)
    }
    ## parallel::stopCluster(cl)

    print("from_list_projected is done")
    names(flow_dist_from_list_projected) <- projected_nodeIDs_vec
    save(flow_dist_from_list_projected,
         file = paste0(file_path,"inference/flow_dist_from_list_projected.RData"))

    load(file = paste0(file_path, "inference/flow_dist_from_list_polluters.RData"))

    # Appending the from list for polluters and projected
    ## flow_dist_from_list_projected = projected_nodeIDs_vec
    flow_dist_from_list <- append(flow_dist_from_list_polluters,flow_dist_from_list_projected)
    # Saving the flow distance from list
    save(flow_dist_from_list,file = paste0(file_path,"inference/flow_dist_from_list.RData"))


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
  }}
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
