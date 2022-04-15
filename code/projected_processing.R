file_path <- "/Users/Amal/Box Sync/PSU/Geoscience_Research/River network project/Analysis for Spill Paper/To_Alex/Cl_spill_whole/revision_analysis/3km/"

missing_ids <- c()
for (i in 1:nrow(df_projected_nodeIDs)){
  trial <- try(load(file = paste0(file_path,"anpoll_files/flow_projected/flow_from_",i,".RData")), silent = TRUE)
  if(!('character' %in% class(trial))){
    missing_ids <- c(missing_ids, i)
  }
}

save(missing_ids, file = paste0(file_path,"anpoll_files/projected_missing_ids.RData"))

##############################################################################################
flow_dist_from_list_projected<-list()
for (i in 1:nrow(df_projected_nodeIDs)){
  load(file = paste0(file_path,"anpoll_files/flow_projected/flow_from_",i,".RData"))
  flow_dist_from_list_projected[[i]] <- flow_dist_vec
  if (i%%1000==0){
    print(i)
  }
}
str(flow_dist_from_list_projected)

names(flow_dist_from_list_projected) <- projected_nodeIDs_vec
print("from_list_projected is done")

save(flow_dist_from_list_projected,file = paste0(file_path,"inference/flow_dist_from_list_projected.RData"))

load(file = paste0(file_path,"inference/flow_dist_from_list_polluters.RData"))
str(flow_dist_from_list_polluters)

flow_dist_from_list <- append(flow_dist_from_list_polluters,flow_dist_from_list_projected)
## Saving the flow distance from list
save(flow_dist_from_list,file = paste0(file_path,"inference/flow_dist_from_list.RData"))





