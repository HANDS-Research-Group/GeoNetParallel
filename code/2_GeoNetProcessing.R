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
df_polluter_raw<-read.csv(file = paste0(file_path,"data/Pollution_Site_R_Package.csv"))
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
df_anpoll_preprocessed %>% as_tibble() %>% group_by(anpoll_indicator) %>% summarise(n = n())

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
    distm(x = as.matrix(df_anpoll_processed[1,c("lon","lat")]),y = as.matrix(df_anpoll_processed[1,c("lon_mapped","lat_mapped")]))

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

cl <- makeCluster(spec = n_chunks)
registerDoParallel(cl)
shortest_path_anpoll_edgelist_chunk <- foreach (index = 1:n_chunks)%dopar% {
    shortest_path_edgelist_creator_parallelized_2(n_chunks=n_chunks,file_path = file_path)
}
stopCluster(cl)

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

                                        # ########################################################################################################
                                        # ########################################################################################################
                                        # ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###
                                        # ## Test code!!
                                        # {## checking vertex ids for shortest path calculations
                                        #   load(file = paste0(file_path,"analyte_files/analyte_vertex_ids.RData"))
                                        #   str(analyte_vertex_ids)
                                        # 
                                        #   load(file = paste0(file_path,"polluter_files/polluter_vertex_ids.RData"))
                                        #   str(polluter_vertex_ids)
                                        # 
                                        #   load(file = paste0(file_path,"polluter_files/df_polluter_processed.RData"))
                                        #   str(df_polluter_processed)
                                        # 
                                        #   load(file = paste0(file_path,"polluter_files/projected_nodeIDs_list.RData"))
                                        #   str(projected_nodeIDs_list)
                                        # 
                                        #   load(file = paste0(file_path,"polluter_files/projected_nodeIDs_vec.RData"))
                                        #   str(projected_nodeIDs_vec)
                                        # 
                                        #   ## correcting the projected_nodeIDs_vec
                                        #   projected_nodeIDs_vec<-unlist(projected_nodeIDs_list)
                                        #   attr(projected_nodeIDs_vec,"names")<-NULL
                                        #   save(projected_nodeIDs_vec,file = paste0(file_path,"polluter_files/projected_nodeIDs_vec.RData"))
                                        # 
                                        #   load(file = paste0(file_path,"polluter_files/projected_nodeIDs_vec.RData"))
                                        #   str(projected_nodeIDs_vec)
                                        # 
                                        #   ## checking the shortest path edgelist
                                        #   load(file = paste0(file_path,"anpoll_files/shortest_path_anpoll_edgelist.RData"))
                                        #   str(shortest_path_anpoll_edgelist)
                                        #   }
                                        # 

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
library(Rcpp)
sourceCpp('flow_dist_cal_cpp.cpp')

## Defining a wrapper function to calculate flow distance of all polluter nodes in parallel
wrapper_flow_dist_cal <- function(df_polluter,anpoll_edgelist, shortest_path_anpoll_edgelist, total_edgelist_character_modified, from_indicator, to_indicator, flow_type, file_path){
    print("index")
    print(index)
    ## Defining a function to calculate between a test node and a
    ## current node. The current nodes are usually polluter nodes.
    flow_dist_cal<-function(anpoll_edgelist,
                            shortest_path_anpoll_edgelist,
                            total_edgelist_character_modified,
                            test_node_ID, current_node_ID,
                            from_indicator, to_indicator, file_path){
        ## Getting the edge row ID
        if(from_indicator){
            edge_row_ID<-which((anpoll_edgelist[,1]==test_node_ID) &
                               (anpoll_edgelist[,2]==current_node_ID))
        }else if(to_indicator){
            edge_row_ID<-which(
            (anpoll_edgelist[,1]==current_node_ID) &
            (anpoll_edgelist[,2]==test_node_ID))
        }
        ## Initializing and loop to fill up the matrix of edges that
        ## are formed by the nodeIDs in the path of test and current
        ## node ID
        flow_path_nodeIDs<- matrix(NA_character_,
                                   nrow =
                                       (length(shortest_path_anpoll_edgelist[[edge_row_ID]])-1),
                                   ncol = 2)
        for (i in 1:(length(shortest_path_anpoll_edgelist[[edge_row_ID]])-1)){
            flow_path_nodeIDs[i,]<- c(shortest_path_anpoll_edgelist[[edge_row_ID]][i],
                                      shortest_path_anpoll_edgelist[[edge_row_ID]][i+1])
        }
        ## Loading the total_edgelist_character_modified
        load(file = paste0(file_path,"common_files_modified/stream_path_dist_vec.RData"))
        ## Initializing the flow distance in meters
        flow_dist_m<-0
        ## Loop to add up the flow distance of all edges in the flow_path_nodeIDs
        ## for(i in 1:nrow(flow_path_nodeIDs)){
        ##     additional_dist <- stream_path_dist_vec[which(
        ##     (total_edgelist_character_modified[,1]==flow_path_nodeIDs[i,1]) &
        ##     (total_edgelist_character_modified[,2]==flow_path_nodeIDs[i,2]))]
        ##     if(length(additional_dist)>0){
        ##         flow_dist_m <- flow_dist_m + as.numeric(stream_path_dist_vec[which(
        ##         (total_edgelist_character_modified[,1]==flow_path_nodeIDs[i,1]) &
        ##         (total_edgelist_character_modified[,2]==flow_path_nodeIDs[i,2]))])
        ##     }
        ## }

        flow_dist_m  <- getDistance(flow_path_nodeIDs = flow_path_nodeIDs,
                                      total_edgelist = total_edgelist_character_modified,
                                      stream_path_dist = stream_path_dist_vec)

        ## flow distance in km
        flow_dist_km<-flow_dist_m/1000
        
        return(flow_dist_km)
    }
    ## Extracting all connected nodeIDs
    if(from_indicator){
        connected_nodeIDs<-anpoll_edgelist[which(anpoll_edgelist[,2]==df_polluter$nodeID[index]),1]
    }else if(to_indicator){
        connected_nodeIDs<-anpoll_edgelist[which(anpoll_edgelist[,1]==df_polluter$nodeID[index]),2]
    }
    ## Initializing the flow_dist_vec as NA. If there are no connected
    ## nodes, this would remain NA
    flow_dist_vec<-NA
    if(length(connected_nodeIDs)>0){
        flow_dist_vec<-rep(NA_integer_,length(connected_nodeIDs))
        print(length(connected_nodeIDs))
        for(j in 1:length(connected_nodeIDs)){
            ## debug(flow_dist_cal)
            flow_dist_vec[j]<-flow_dist_cal(anpoll_edgelist = anpoll_edgelist,
                                            shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist,
                                            total_edgelist_character_modified = total_edgelist_character_modified,
                                            test_node_ID = connected_nodeIDs[j],
                                            current_node_ID = df_polluter$nodeID[index],
                                            from_indicator = from_indicator,
                                            to_indicator = to_indicator,
                                            file_path = file_path)
            ## print(j)
        }
    }
    save(flow_dist_vec, file = paste0(file_path,"anpoll_files/",
                                      flow_type, "/flow_from_", index,
                                      ".RData"))
    ## save(flow_dist_vec, file = paste0(file_path,"anpoll_files/projected_missing_ids/flow_from_", index,".RData"))
    return(flow_dist_vec)
}

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
     cl <- makeCluster(39)
     registerDoParallel(cl)
     flow_dist_from_list<-foreach (index = 1:nrow(df_polluter_nodeID_aggregated_rank_subgaph))%dopar% {
         wrapper_flow_dist_cal(df_polluter=df_polluter_nodeID_aggregated_rank_subgaph,anpoll_edgelist = anpoll_edgelist, shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist, total_edgelist_character_modified = total_edgelist_character_modified, from_indicator = T, to_indicator = F, file_path = file_path)
     }
     stopCluster(cl)
     names(flow_dist_from_list)<-df_polluter_nodeID_aggregated_rank_subgaph$nodeID
     save(flow_dist_from_list,file = paste0(file_path,"inference/flow_dist_from_list.RData"))
     print("from_list is done")
     ## Getting the flow dist to list
     flow_dist_to_list<-list()
     cl <- makeCluster(39)
     registerDoParallel(cl)
     flow_dist_to_list<-foreach (index = 1:nrow(df_polluter_nodeID_aggregated_rank_subgaph))%dopar% {
         wrapper_flow_dist_cal(df_polluter=df_polluter_nodeID_aggregated_rank_subgaph,anpoll_edgelist = anpoll_edgelist, shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist, total_edgelist_character_modified = total_edgelist_character_modified, from_indicator = F, to_indicator = T, file_path = file_path)
     }
     stopCluster(cl)
     names(flow_dist_to_list)<-df_polluter_nodeID_aggregated_rank_subgaph$nodeID
     save(flow_dist_to_list,file = paste0(file_path,"inference/flow_dist_to_list.RData"))
 }else{
     ## Getting the flow dist from list for polluters
     flow_dist_from_list_polluters<-list()
     ## for (index in 1:nrow(df_polluter_nodeID_aggregated)){
     ##     flow_dist_from_list_polluters[[index]] <- wrapper_flow_dist_cal(df_polluter=df_polluter_nodeID_aggregated,anpoll_edgelist = anpoll_edgelist, shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist, total_edgelist_character_modified = total_edgelist_character_modified, from_indicator = T, to_indicator = F, file_path = file_path)
     ## }
     ## index <- 1401
     ## debug(wrapper_flow_dist_cal)
     ## temp <- wrapper_flow_dist_cal(df_polluter=df_polluter_nodeID_aggregated,anpoll_edgelist = anpoll_edgelist, shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist, total_edgelist_character_modified = total_edgelist_character_modified, from_indicator = T, to_indicator = F, flow_type = "flow_from", file_path = file_path)

     ## ptm <- proc.time()
     print("From_list_polluter started")

     ## Uncomment the following section
     # cl <- makeCluster(30)
     # registerDoParallel(cl)
     ## flow_dist_from_list_polluters <- foreach (index = 1:nrow(df_polluter_nodeID_aggregated)) %dopar%{
     flow_dist_from_list_polluters <- foreach (index = 1:nrow(df_polluter_nodeID_aggregated), .packages="Rcpp", .noexport = "getDistance") %do% {

         sourceCpp('flow_dist_cal_cpp.cpp')
         wrapper_flow_dist_cal(df_polluter = df_polluter_nodeID_aggregated,
                               anpoll_edgelist = anpoll_edgelist,
                               shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist,
                               total_edgelist_character_modified =
                                   total_edgelist_character_modified,
                               from_indicator = T, to_indicator = F, flow_type = "flow_from",
                               file_path = file_path)
     }
     # stopCluster(cl)
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

     ## flow_dist_from_list_projected  <- vector(mode="list", length = length(projected_nodeIDs_vec))
     ## index <- 4
     ## flow_dist_from_list_projected[1]  <- wrapper_flow_dist_cal(df_polluter=df_projected_nodeIDs,
     ##                           anpoll_edgelist = anpoll_edgelist,
     ##                           shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist,
     ##                           total_edgelist_character_modified = total_edgelist_character_modified,
     ##                           from_indicator = T, to_indicator = F,
     ##                           flow_type = "flow_projected",
     ##                           file_path = file_path)

     # cl <- makeCluster(50)
     # registerDoParallel(cl)

     flow_dist_from_list_projected<-foreach (index =1:nrow(df_projected_nodeIDs), .packages="Rcpp", .noexport = "getDistance")%do% {
         sourceCpp('flow_dist_cal_cpp.cpp')
         wrapper_flow_dist_cal(df_polluter=df_projected_nodeIDs,
                               anpoll_edgelist = anpoll_edgelist,
                               shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist,
                               total_edgelist_character_modified = total_edgelist_character_modified,
                               from_indicator = T, to_indicator = F,
                               flow_type = "flow_projected",
                               file_path = file_path)
     }
     # stopCluster(cl)

     print("from_list_projected is done")
     print(length(projected_nodeIDs_vec))
     print(length(flow_dist_from_list_projected))
     names(flow_dist_from_list_projected) <- projected_nodeIDs_vec
     save(flow_dist_from_list_projected,
          file = paste0(file_path,"inference/flow_dist_from_list_projected.RData"))

     ## ## Initializing the combined list
     ## flow_dist_from_list_projected_combined<-list()
     ## ## for (i in 1:100){
     ## ##     load(file = paste0(file_path,"inference/projected_results/flow_dist_from_list_projected_",i,".RData"))
     ## ##     flow_dist_from_list_projected_combined[[i]]<-flow_dist_from_list_projected
     ## ## }
     
     ##                                    #str(flow_dist_from_list_projected_combined[[100]])
     ## flow_dist_from_list_projected<-unlist(x = flow_dist_from_list_projected_combined,
     ##                                       recursive = F)
     ## str(flow_dist_from_list_projected)


     load(file = paste0(file_path, "inference/flow_dist_from_list_polluters.RData"))

     # Appending the from list for polluters and projected
     ## flow_dist_from_list_projected = projected_nodeIDs_vec
     flow_dist_from_list <- append(flow_dist_from_list_polluters,flow_dist_from_list_projected)
     # Saving the flow distance from list
     save(flow_dist_from_list,file = paste0(file_path,"inference/flow_dist_from_list.RData"))


     ## Getting the flow dist to list
     flow_dist_to_list<-list()

     flow_dist_to_list  <- vector(mode="list", length=length( df_polluter_nodeID_aggregated))
     ## flow_dist_to_list[1] <- wrapper_flow_dist_cal(df_polluter= df_polluter_nodeID_aggregated,
     ##                           anpoll_edgelist = anpoll_edgelist,
     ##                           shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist,
     ##                           total_edgelist_character_modified =
     ##                               total_edgelist_character_modified,
     ##                           from_indicator = F, to_indicator = T,
     ##                           flow_type = "flow_to", file_path = file_path)

     ## Actual Code
     # cl <- makeCluster(30)
     # registerDoParallel(cl)
     flow_dist_to_list<-foreach (index = 1:nrow(df_polluter_nodeID_aggregated), .packages="Rcpp", .noexport = "getDistance")%do% {
         sourceCpp('flow_dist_cal_cpp.cpp')
         wrapper_flow_dist_cal(df_polluter= df_polluter_nodeID_aggregated,
                               anpoll_edgelist = anpoll_edgelist,
                               shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist,
                               total_edgelist_character_modified =
                                   total_edgelist_character_modified,
                               from_indicator = F, to_indicator = T,
                               flow_type = "flow_to", file_path = file_path)
     }
     # stopCluster(cl)

     names(flow_dist_to_list)<-df_polluter_nodeID_aggregated$nodeID
     print(" flow_dist_to_list is done")
     save(flow_dist_to_list,file = paste0(file_path,"inference/flow_dist_to_list.RData"))
 }}

