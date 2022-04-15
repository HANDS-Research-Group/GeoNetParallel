file_path <- "/Users/Amal/Box Sync/PSU/Geoscience_Research/River network project/Analysis for Spill Paper/To_Alex/Cl_spill_whole/revision_analysis/"

########################################################################################################
########################################################################################################
########################################################################################################
############################################# Inference ################################################
########################################################################################################
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
wrapper_flow_dist_cal <- function(df_polluter,anpoll_edgelist, shortest_path_anpoll_edgelist, total_edgelist_character_modified, from_indicator, to_indicator, file_path){
  ## Defining a function to calculate between a test node and a current node. The current nodes are usually polluter nodes.
  flow_dist_cal<-function(anpoll_edgelist, shortest_path_anpoll_edgelist, total_edgelist_character_modified, test_node_ID, current_node_ID, from_indicator, to_indicator, file_path){
    ## Getting the edge row ID
    if(from_indicator){
      edge_row_ID<-which((anpoll_edgelist[,1]==test_node_ID)&(anpoll_edgelist[,2]==current_node_ID))
    }else if(to_indicator){
      edge_row_ID<-which((anpoll_edgelist[,1]==current_node_ID)&(anpoll_edgelist[,2]==test_node_ID))
    }
    ## Initializing and loop to fill up the matrix of edges that are formed by the nodeIDs in the path of test and current node ID
    flow_path_nodeIDs<-matrix(NA_character_,nrow = (length(shortest_path_anpoll_edgelist[[edge_row_ID]])-1),ncol = 2)
    for (i in 1:(length(shortest_path_anpoll_edgelist[[edge_row_ID]])-1)){
      flow_path_nodeIDs[i,]<-c(shortest_path_anpoll_edgelist[[edge_row_ID]][i],shortest_path_anpoll_edgelist[[edge_row_ID]][i+1])
    }
    ## Loading the total_edgelist_character_modified
    load(file = paste0(file_path,"common_files_modified/stream_path_dist_vec.RData"))
    ## Initializing the flow distance in meters
    flow_dist_m<-0
    ## Loop to add up the flow distance of all edges in the flow_path_nodeIDs
    for(i in 1:nrow(flow_path_nodeIDs)){
      additional_dist <- stream_path_dist_vec[which((total_edgelist_character_modified[,1]==flow_path_nodeIDs[i,1])&(total_edgelist_character_modified[,2]==flow_path_nodeIDs[i,2]))]
      if(length(additional_dist)>0){
        flow_dist_m <- flow_dist_m + as.numeric(stream_path_dist_vec[which((total_edgelist_character_modified[,1]==flow_path_nodeIDs[i,1])&(total_edgelist_character_modified[,2]==flow_path_nodeIDs[i,2]))])
      }
    }
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
  ## Initializing the flow_dist_vec as NA. If there are no connected nodes, this would remain NA
  flow_dist_vec<-NA
  if(length(connected_nodeIDs)>0){
    flow_dist_vec<-rep(NA_integer_,length(connected_nodeIDs))
    for(j in 1:length(connected_nodeIDs)){
      #debug(flow_dist_cal)
      flow_dist_vec[j]<-flow_dist_cal(anpoll_edgelist = anpoll_edgelist,shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist,total_edgelist_character_modified = total_edgelist_character_modified,test_node_ID = connected_nodeIDs[j],current_node_ID = df_polluter$nodeID[index],from_indicator = from_indicator,to_indicator = to_indicator,file_path = file_path)
      print(j)
    }
  }
  return(flow_dist_vec)
}

## Loading the dataframe "df_anpoll_processed"
load(file = paste0(file_path, "anpoll_files/df_anpoll_processed.RData"))
## Loading the vector of projected node IDs
load(file = paste0(file_path, "polluter_files/projected_nodeIDs_vec.RData"))


## Subsetting the df_polluter_processed from df_anpoll_processed
df_polluter_processed <- df_anpoll_processed[which(df_anpoll_processed$anpoll_indicator == "polluter"),]
nrow(df_polluter_processed)

flow_dist_from_list_polluters <- list()

for (index in 1:nrow(df_polluter_nodeID_aggregated)){
  flow_dist_from_list_polluters[[index]] <- wrapper_flow_dist_cal(df_polluter=df_polluter_nodeID_aggregated,anpoll_edgelist = anpoll_edgelist, shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist, total_edgelist_character_modified = total_edgelist_character_modified, from_indicator = T, to_indicator = F, file_path = file_path)
  print(index)
}


debug(wrapper_flow_dist_cal)
index <- 47
flow_dist_from_list_polluters[[index]] <- wrapper_flow_dist_cal(df_polluter=df_polluter_nodeID_aggregated,anpoll_edgelist = anpoll_edgelist, shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist, total_edgelist_character_modified = total_edgelist_character_modified, from_indicator = T, to_indicator = F, file_path = file_path)



flow_dist_from_list_polluters<-list()
cl <- makeCluster(19)
registerDoParallel(cl)
flow_dist_from_list_polluters<-foreach (index = 1:nrow(df_polluter_nodeID_aggregated))%dopar% {
  wrapper_flow_dist_cal(df_polluter=df_polluter_nodeID_aggregated,anpoll_edgelist = anpoll_edgelist, shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist, total_edgelist_character_modified = total_edgelist_character_modified, from_indicator = T, to_indicator = F, file_path = file_path)
}
stopCluster(cl)
names(flow_dist_from_list_polluters)<-df_polluter_nodeID_aggregated$nodeID
save(flow_dist_from_list_polluters,file = paste0(file_path,"inference/flow_dist_from_list_polluters.RData"))
print("from_list_polluter is done")