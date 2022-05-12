###############################################################################################
################################## Cl Spill Analysis ##########################################
###############################################################################################
###############################################################################################

library(tidyverse)

## Update the path to the folder
file_path <- "C:/GeoNet/GeoNet_2021_packageDataset/"

## Sourcing the modular functions for analysis
source(file = paste0(file_path, "code/BaSu_network_v16_modular_functions.R"))

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
specCond <- read.csv(file = paste0(file_path, 'data/Water_Chemistry_R_Package.csv'))
df_analyte_preprocessed = specCond[,c(1, 2, 3, 4)]
names(df_analyte_preprocessed) <- c("conc", "lon", "lat", "date")

df_analyte_preprocessed$date <- as.Date(df_analyte_preprocessed$date, format = "%m/%d/%y")


## analyte_data <- read.csv(file = paste0(file_path, 'data/Cl_data.csv'))
## df_analyte_preprocessed = analyte_data[,c(1, 2, 3)]
## df_analyte_preprocessed$date = as.Date(analyte_data$Date, "%m/%d/%y")
## names(df_analyte_preprocessed) <- c("conc", "lon", "lat", "date")
str(df_analyte_preprocessed)
head(df_analyte_preprocessed)
## save(df_analyte_preprocessed, file=paste0(file_path, "analyte_files/df_analyte_preprocessed.RData"))

###################################################
## Preprocessing the polluter dataframe
df_polluter_raw<-read.csv(file = paste0(file_path,"data/Pollution_Site_R_Package_missing_1.csv"))
str(df_polluter_raw)

## Preprocessing the analyte dataframe
df_polluter_preprocessed<-df_polluter_raw[,c(2,3)]
df_polluter_preprocessed$date<-as.Date(df_polluter_raw$date,format="%m/%d/%Y")
str(df_polluter_preprocessed)
head(df_polluter_preprocessed)
## saving this dataframe
save(df_polluter_preprocessed,file = paste0(file_path,"polluter_files/df_polluter_preprocessed.RData"))

####################################################
## Combining the analyte and polluter dataframes
df_anpoll_preprocessed <- data.frame(matrix(NA,nrow(df_analyte_preprocessed)+nrow(df_polluter_preprocessed),4))
df_anpoll_preprocessed <- data.frame(df_analyte_preprocessed,"anpoll_indicator"=rep("analyte",nrow(df_analyte_preprocessed)),stringsAsFactors = F)
df_anpoll_preprocessed[(nrow(df_analyte_preprocessed)+1):(nrow(df_analyte_preprocessed)+nrow(df_polluter_preprocessed)),c(2,3,4)] <- df_polluter_preprocessed
df_anpoll_preprocessed[(nrow(df_analyte_preprocessed)+1):(nrow(df_analyte_preprocessed)+nrow(df_polluter_preprocessed)),"anpoll_indicator"] <- "polluter"

str(df_anpoll_preprocessed)
head(df_anpoll_preprocessed)
tail(df_anpoll_preprocessed)

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
{## Test code!!
    ## checking the output files of two step anpoll mapper function and others
    load(file=paste0(file_path,"anpoll_files/dist_mapped.RData"))
    str(dist_mapped)
    summary(dist_mapped)
    hist(dist_mapped)
    hist(dist_mapped[dist_mapped<200])

    load(file=paste0(file_path,"anpoll_files/df_anpoll_processed.RData"))
    str(df_anpoll_processed)
    head(df_anpoll_processed)
    length(unique(df_anpoll_processed$nodeID))
    geosphere::distm(x = as.matrix(df_anpoll_processed[1,c("lon","lat")]),y = as.matrix(df_anpoll_processed[1,c("lon_mapped","lat_mapped")]))

    load(file=paste0(file_path,"polluter_files/df_polluter_processed.RData"))
    str(df_polluter_processed)
    summary(tail(dist_mapped))

    load(file=paste0(file_path,"common_files_modified/df_node_latlong_modified.RData"))
    str(df_node_latlong_modified)
    head(df_node_latlong_modified)

    load(file=paste0(file_path,"common_files/df_node_latlong.RData"))
    str(df_node_latlong)
    length(unique(df_node_latlong$nodeID))

    load(file = paste0(file_path,"common_files/stream_list.RData"))
    str(stream_list)

    load(file = paste0(file_path,"common_files_modified/stream_list_modified.RData"))
    str(stream_list_modified)

    load(file = paste0(file_path,"common_files_modified/stream_path_dist_vec.RData"))
    str(stream_path_dist_vec)
    summary(stream_path_dist_vec)

    load(file = paste0(file_path,"analyte_files/df_analyte_nodeID_aggregated.RData"))
    str(df_analyte_nodeID_aggregated)

    load(file = paste0(file_path,"analyte_files/list_analyte_time.RData"))
    str(list_analyte_time, max.level = 1)

    load(file = paste0(file_path,"analyte_files/listID_nodeID_matrix.RData"))
    str(listID_nodeID_matrix)
    head(listID_nodeID_matrix)

    load(file = paste0(file_path,"polluter_files/df_polluter_nodeID_aggregated.RData"))
    str(df_polluter_nodeID_aggregated)

    load(file = paste0(file_path,"common_files_modified/igraph_river_decomposed_list.RData"))
    str(igraph_river_decomposed_list, max.level = 1)}

##  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###
#######################################################################################################

## Getting the analyte_vertex_ids
analyte_vertex_ids <- vertex_IDs_generator(df_anpoll_nodeID_aggregated = df_analyte_nodeID_aggregated,output_path = file_path,analyte=T,polluter=F,graph_subset=graph_subset_indicator)

## Getting the polluter_vertex_ids
polluter_vertex_ids <- vertex_IDs_generator(df_anpoll_nodeID_aggregated = df_polluter_nodeID_aggregated,output_path = file_path,analyte=F,polluter=T, graph_subset=graph_subset_indicator)

## Getting the projected_nodeIDs_list
projected_nodeIDs_list <- projected_nodeIDs_list_generator(file_path=file_path, projected_threshold_dist_km = 50)

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

## Defining number of chunks to be parallelized, chunk size and indices inside each chunk
n_chunks <- 4
chunk_size <- ceiling(length(analyte_polluter_projected_vertex_ids)/n_chunks)
indices_all <- 1:length(analyte_polluter_projected_vertex_ids)
indices_chunk_list <- split(x = indices_all, ceiling(seq_along(indices_all)/chunk_size))

cl <- parallel::makeCluster(spec = n_chunks)
doParallel::registerDoParallel(cl)
shortest_path_anpoll_edgelist_chunk <- foreach (index = 1:n_chunks)%dopar% {
    shortest_path_edgelist_creator_parallelized_2(n_chunks=n_chunks,file_path = file_path)
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

head(anpoll_network_list[[1]])

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
     cl <- parallel::makeCluster(39)
     doParallel::registerDoParallel(cl)
     flow_dist_from_list<-foreach (index = 1:nrow(df_polluter_nodeID_aggregated_rank_subgaph))%dopar% {
         wrapper_flow_dist_cal(index=index, df_polluter=df_polluter_nodeID_aggregated_rank_subgaph,anpoll_edgelist = anpoll_edgelist, shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist, total_edgelist_character_modified = total_edgelist_character_modified, from_indicator = T, to_indicator = F, file_path = file_path)
     }
     parallel::stopCluster(cl)
     names(flow_dist_from_list)<-df_polluter_nodeID_aggregated_rank_subgaph$nodeID
     save(flow_dist_from_list,file = paste0(file_path,"inference/flow_dist_from_list.RData"))
     print("from_list is done")
     ## Getting the flow dist to list
     flow_dist_to_list<-list()
     cl <- parallel::makeCluster(39)
     doParallel::registerDoParallel(cl)
     flow_dist_to_list<-foreach (index = 1:nrow(df_polluter_nodeID_aggregated_rank_subgaph))%dopar% {
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
     # cl <- parallel::makeCluster(30)
     # doParallel::registerDoParallel(cl)
     ## flow_dist_from_list_polluters <- foreach (index = 1:nrow(df_polluter_nodeID_aggregated)) %dopar%{
     sourceCpp( paste0(file_path, "code/flow_dist_cal_cpp.cpp") )
     flow_dist_from_list_polluters <- foreach (index = 1:nrow(df_polluter_nodeID_aggregated), .packages="Rcpp", .noexport = "getDistance") %do% {

         #sourceCpp( paste0(file_path, "code/flow_dist_cal_cpp.cpp") )
         wrapper_flow_dist_cal(index, df_polluter = df_polluter_nodeID_aggregated,
                               anpoll_edgelist = anpoll_edgelist,
                               shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist,
                               total_edgelist_character_modified =
                                   total_edgelist_character_modified,
                               from_indicator = T, to_indicator = F, flow_type = "flow_from",
                               file_path = file_path)
     }
     # parallel::stopCluster(cl)
     print(length(flow_dist_from_list_polluters))
     print(length(df_polluter_nodeID_aggregated$nodeID))
     names(flow_dist_from_list_polluters)<-df_polluter_nodeID_aggregated$nodeID
     save(flow_dist_from_list_polluters,
          file = paste0(file_path,"inference/flow_dist_from_list_polluters.RData"))
     print("from_list_polluter is done")


     projected_nodeIDs_vec = unique(projected_nodeIDs_vec)
     ## Getting the flow dist from list for projected node IDs of polluters
     df_projected_nodeIDs<-data.frame("nodeID"=projected_nodeIDs_vec,stringsAsFactors = F)
     flow_dist_from_list_projected<-list()

     # cl <- parallel::makeCluster(50)
     # doParallel::registerDoParallel(cl)
    sourceCpp( paste0(file_path, "code/flow_dist_cal_cpp.cpp") )
     flow_dist_from_list_projected<-foreach (index =1:nrow(df_projected_nodeIDs), .packages="Rcpp", .noexport = "getDistance")%do% {
         sourceCpp('flow_dist_cal_cpp.cpp')
         wrapper_flow_dist_cal(index, df_polluter=df_projected_nodeIDs,
                               anpoll_edgelist = anpoll_edgelist,
                               shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist,
                               total_edgelist_character_modified = total_edgelist_character_modified,
                               from_indicator = T, to_indicator = F,
                               flow_type = "flow_projected",
                               file_path = file_path)
     }
     # parallel::stopCluster(cl)

     print("from_list_projected is done")
     print(length(projected_nodeIDs_vec))
     print(length(flow_dist_from_list_projected))
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
     # cl <- parallel::makeCluster(30)
     # doParallel::registerDoParallel(cl)
     sourceCpp( paste0(file_path, "code/flow_dist_cal_cpp.cpp") )
     flow_dist_to_list<-foreach (index = 1:nrow(df_polluter_nodeID_aggregated), .packages="Rcpp", .noexport = "getDistance")%do% {
         sourceCpp('flow_dist_cal_cpp.cpp')
         wrapper_flow_dist_cal(index, df_polluter= df_polluter_nodeID_aggregated,
                               anpoll_edgelist = anpoll_edgelist,
                               shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist,
                               total_edgelist_character_modified =
                                   total_edgelist_character_modified,
                               from_indicator = F, to_indicator = T,
                               flow_type = "flow_to", file_path = file_path)
     }
     # parallel::stopCluster(cl)

     names(flow_dist_to_list)<-df_polluter_nodeID_aggregated$nodeID
     print(" flow_dist_to_list is done")
     save(flow_dist_to_list,file = paste0(file_path,"inference/flow_dist_to_list.RData"))
 }}

