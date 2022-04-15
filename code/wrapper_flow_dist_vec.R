library(tidyverse)
library(Rcpp)
## Update the path to the folder
file_path <- "/home/rrpatil/GeoNet_2021_test/"
# file_path <- "/gpfs/group/eesi/default/users/aua257/revision_analysis/3km/"
# file_path <- "/Users/Amal/Box Sync/PSU/Geoscience_Research/River network project/Analysis for Spill Paper/To_Alex/Cl_spill_whole/revision_analysis/3km/"

## Sourcing the modular functions for analysis
source(file = paste0(file_path, "code/BaSu_network_v16_modular_functions.R"))
sourceCpp('flow_dist_cal_cpp.cpp')
################################################################################
########################################### Inference ##########################
################################################################################
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
{
    if(graph_subset_indicator)
    {
        df_polluter_nodeID_aggregated_rank_subgaph<-df_polluter_nodeID_aggregated[which(df_polluter_nodeID_aggregated$rank_subgraph),]
        row.names(df_polluter_nodeID_aggregated_rank_subgaph)<-NULL
        ## Getting the flow dist from list
        flow_dist_from_list<-list()
        cl <- makeCluster(10)
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
        cl <- makeCluster(10)
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


     ## Uncomment the following section
     cl <- makeCluster(50)
     registerDoParallel(cl)
     ## flow_dist_from_list_polluters <- foreach (index = 1:nrow(df_polluter_nodeID_aggregated)) %dopar%{
     flow_dist_from_list_polluters <- foreach (index = 1:nrow(df_polluter_nodeID_aggregated), .packages="Rcpp", .noexport = "getDistance") %dopar% {

         sourceCpp('flow_dist_cal_cpp.cpp')
         wrapper_flow_dist_cal(df_polluter = df_polluter_nodeID_aggregated,
                               anpoll_edgelist = anpoll_edgelist,
                               shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist,
                               total_edgelist_character_modified =
                                   total_edgelist_character_modified,
                               from_indicator = T, to_indicator = F, flow_type = "flow_from",
                               file_path = file_path)
     }
     stopCluster(cl)
     print(length(flow_dist_from_list_polluters))
     print(length(df_polluter_nodeID_aggregated$nodeID))
     names(flow_dist_from_list_polluters)<-df_polluter_nodeID_aggregated$nodeID
     save(flow_dist_from_list_polluters,
          file = paste0(file_path,"inference/flow_dist_from_list_polluters.RData"))
     print("from_list_polluter is done")
     ## print("Time : ")
     ## print(proc.time()-ptm)


        ## Getting the flow dist from list for projected node IDs of polluters
        print(length(projected_nodeIDs_vec))
        projected_nodeIDs_vec = unique(projected_nodeIDs_vec)
        print(length(projected_nodeIDs_vec))
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

     cl <- makeCluster(40)
     registerDoParallel(cl)

     ## flow_dist_from_list_projected<-foreach (index = 1:10)%dopar% {
     flow_dist_from_list_projected<-foreach (index =1:nrow(df_projected_nodeIDs), .packages="Rcpp", .noexport = "getDistance")%dopar% {
         sourceCpp('flow_dist_cal_cpp.cpp')
         wrapper_flow_dist_cal(df_polluter=df_projected_nodeIDs,
                               anpoll_edgelist = anpoll_edgelist,
                               shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist,
                               total_edgelist_character_modified = total_edgelist_character_modified,
                               from_indicator = T, to_indicator = F,
                               flow_type = "flow_projected",
                               file_path = file_path)
     }
     stopCluster(cl)

     print("from_list_projected is done")
     print(length(projected_nodeIDs_vec))
     print(length(flow_dist_from_list_projected))
     names(flow_dist_from_list_projected) <- projected_nodeIDs_vec
     save(flow_dist_from_list_projected,
          file = paste0(file_path,"inference/flow_dist_from_list_projected.RData"))
     ## print("Time : ")
     ## print(proc.time()-ptm)

     ## Initializing the combined list
     flow_dist_from_list_projected_combined<-list()
     ## for (i in 1:100){
     ##     load(file = paste0(file_path,"inference/projected_results/flow_dist_from_list_projected_",i,".RData"))
     ##     flow_dist_from_list_projected_combined[[i]]<-flow_dist_from_list_projected
     ## }
     
                                        #str(flow_dist_from_list_projected_combined[[100]])
     ## flow_dist_from_list_projected<-unlist(x = flow_dist_from_list_projected_combined,
     ##                                       recursive = F)
     ## str(flow_dist_from_list_projected)


     ## load(file = paste0(file_path, "inference/flow_dist_from_list_polluters.RData"))

     # Appending the from list for polluters and projected
     names(flow_dist_from_list_projected) <- projected_nodeIDs_vec
     flow_dist_from_list <- append(flow_dist_from_list_polluters,flow_dist_from_list_projected)
     # Saving the flow distance from list
     save(flow_dist_from_list,file = paste0(file_path,"inference/flow_dist_from_list.RData"))

     ## ptm <- proc.time()
     ## Getting the flow dist to list
     flow_dist_to_list<-list()

     flow_dist_to_list  <- vector(mode="list", length=length( df_polluter_nodeID_aggregated))

     ## Actual Code
     cl <- makeCluster(40)
     registerDoParallel(cl)
     flow_dist_to_list<-foreach (index = 1:nrow(df_polluter_nodeID_aggregated), .packages="Rcpp", .noexport = "getDistance")%dopar% {
         sourceCpp('flow_dist_cal_cpp.cpp')
         wrapper_flow_dist_cal(df_polluter= df_polluter_nodeID_aggregated,
                               anpoll_edgelist = anpoll_edgelist,
                               shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist,
                               total_edgelist_character_modified =
                                   total_edgelist_character_modified,
                               from_indicator = F, to_indicator = T,
                               flow_type = "flow_to", file_path = file_path)
     }
     stopCluster(cl)

     names(flow_dist_to_list)<-df_polluter_nodeID_aggregated$nodeID
     print(" flow_dist_to_list is done")
     save(flow_dist_to_list,file = paste0(file_path,"inference/flow_dist_to_list.RData"))
     ## print("Time : ")
     ## print(proc.time()-ptm)

     ## load(file = paste0(file_path,"anpoll_files/projected_missing_ids.RData"))
     ## df_projected_nodeIDs<-data.frame("nodeID"=projected_nodeIDs_vec,stringsAsFactors = F)
     ## ## Getting the flow dist from list for projected node IDs of polluters
     ## cl <- makeCluster(19)
     ## registerDoParallel(cl)
     ## flow_dist_from_list_projected<-foreach (index = missing_ids)%dopar% {
     ##     wrapper_flow_dist_cal(df_polluter=df_projected_nodeIDs,anpoll_edgelist = anpoll_edgelist, shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist, total_edgelist_character_modified = total_edgelist_character_modified, from_indicator = T, to_indicator = F, flow_type = "flow_projected", file_path = file_path)
     ## }
     ## stopCluster(cl)
     ## names(flow_dist_from_list_projected) <- projected_nodeIDs_vec
     ## print("from_list_projected is done")
     
 }}
