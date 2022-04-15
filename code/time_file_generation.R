
library(tidyverse)

## Update the path to the folder
file_path <- "/home/rrpatil/GeoNet_2021_pro/"

## file_path <- "/gpfs/group/eesi/default/users/aua257/revision_analysis/3km/"

## file_path <- "/Users/Amal/Box Sync/PSU/Geoscience_Research/River network project/Analysis for Spill Paper/To_Alex/Cl_spill_whole/revision_analysis/3km/"

## Sourcing the modular functions for analysis
source(file = paste0(file_path, "code/BaSu_network_v16_modular_functions.R"))

## dir.create(path = paste0(file_path,"common_files"))
## dir.create(path = paste0(file_path,"common_files_modified"))
## ## dir.create(path = paste0(file_path,"analyte_files"))
## dir.create(path = paste0(file_path,"polluter_files"))
## dir.create(path = paste0(file_path,"anpoll_files"))
## dir.create(path = paste0(file_path,"anpoll_files/flow_from"))
## dir.create(path = paste0(file_path,"anpoll_files/flow_projected"))
## dir.create(path = paste0(file_path,"anpoll_files/flow_to"))
## dir.create(path = paste0(file_path,"anpoll_files/shortest_path_chunks"))
## dir.create(path = paste0(file_path,"anpoll_files/projected_missing_ids"))
## dir.create(path = paste0(file_path,"inference"))
## dir.create(path = paste0(file_path,"inference/test_results_j1"))
## dir.create(path = paste0(file_path,"inference/test_results_j2"))



###############################################################################################
## Loading the preprocessed analyte dataframe
load(file = paste0(file_path,"analyte_files/df_analyte_preprocessed.RData"))
min(df_analyte_preprocessed$date)
max(df_analyte_preprocessed$date)
## str(df_analyte_preprocessed)


###################################################
## Preprocessing the polluter dataframe
## df_polluter_raw<-read.csv(file = paste0(file_path,"data/newPolluters.csv"),stringsAsFactors = F)
## str(df_polluter_raw)

## ## Preprocessing the analyte dataframe
## df_polluter_preprocessed<-df_polluter_raw[,c(3,2)]
## df_polluter_preprocessed$date<-as.Date(df_polluter_raw$date,format="%m/%d/%Y")
## str(df_polluter_preprocessed)
load(file = paste0(file_path,"polluter_files/df_polluter_preprocessed.RData"))
head(df_polluter_preprocessed)
min(df_polluter_preprocessed$date)
max(df_polluter_preprocessed$date)

## ## saving this dataframe
## save(df_polluter_preprocessed,file = paste0(file_path,"polluter_files/df_polluter_preprocessed.RData"))

####################################################

## Combining the analyte and polluter dataframes
df_anpoll_preprocessed <- data.frame(matrix(NA,nrow(df_analyte_preprocessed)+nrow(df_polluter_preprocessed),4))
df_anpoll_preprocessed <- data.frame(df_analyte_preprocessed,"anpoll_indicator"=rep("analyte",nrow(df_analyte_preprocessed)),stringsAsFactors = F)
df_anpoll_preprocessed[(nrow(df_analyte_preprocessed)+1):(nrow(df_analyte_preprocessed)+nrow(df_polluter_preprocessed)),c(2,3,4)] <- df_polluter_preprocessed
df_anpoll_preprocessed[(nrow(df_analyte_preprocessed)+1):(nrow(df_analyte_preprocessed)+nrow(df_polluter_preprocessed)),"anpoll_indicator"] <- "polluter"

str(df_anpoll_preprocessed)
head(df_anpoll_preprocessed)
tail(df_anpoll_preprocessed)

load(file = paste0(file_path,"anpoll_files/df_anpoll_processed.RData"))

df_anpoll_processed$date  <- df_anpoll_preprocessed$date

save(df_anpoll_processed, file = paste0(file_path,"anpoll_files/df_anpoll_processed.RData"))
df_polluter_processed<-df_anpoll_processed[which(df_anpoll_processed$anpoll_indicator=="polluter"),]
row.names(df_polluter_processed)<-NULL
nrow(df_polluter_processed)
save(df_polluter_processed,file=paste0(file_path,"polluter_files/df_polluter_processed.RData"))

## Getting the 2nd MOST IMPORTANT dataframe and other lists for analyte
df_list_analyte_obj <- df_list_analyte_generator(df_analyte_processed = df_anpoll_processed[which(df_anpoll_processed$anpoll_indicator=="analyte"),],output_path = file_path)

## Getting the 3rd MOST IMPORTANT dataframe for polluter
df_polluter_nodeID_aggregated <- df_polluter_generator(df_polluter_processed = df_anpoll_processed[which(df_anpoll_processed$anpoll_indicator=="polluter"), c("nodeID","lon_mapped","lat_mapped","date")], output_path = file_path)
