load('../polluter_files/projected_nodeIDs_vec.RData')
load('../anpoll_files/anpoll_edgelist.RData')
df_projected_nodeIDs<-data.frame("nodeID"=projected_nodeIDs_vec,stringsAsFactors = F)

library(foreach)
library(doParallel)

cl <- makeCluster(40)
registerDoParallel(cl)

connected_nodes <- foreach(index = 1:nrow(df_projected_nodeIDs)) %dopar% {
  anpoll_edgelist[which(anpoll_edgelist[,2]==df_projected_nodeIDs$nodeID[index]), 2]
}
stopCluster(cl)
str(connected_nodes)

save(connected_nodes, file = 'connected_nodes.RData')
