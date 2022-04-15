library(tidyverse)
print("At the start")
## Update the path to the folder
file_path <- "/home/rrpatil/GeoNet_2021_test/"

## Sourcing the modular functions for analysis
source(file = paste0(file_path, "code/BaSu_network_v16_modular_functions.R"))


############################################# Inference #########################################################################################################################
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


# ########################################################################################################
# ## Summarizing the overall downstream distances for each polluter
# load(file = paste0(file_path,"inference/flow_dist_to_list.RData"))
# flow_dist_to_summary_list <- lapply(X = flow_dist_to_list, FUN = function(x){if (!is.na(x)) return (summary(x)) else return(NULL)})
# flow_dist_to_summary_list[sapply(flow_dist_to_summary_list, is.null)] <- NULL
# df_to_summary <- do.call(rbind,flow_dist_to_summary_list)
# colMeans(df_to_summary)
# mean(df_to_summary[,"Median"])
# 
# ########################################################################################################
# ########################################################################################################
# ## Getting the final matrices and dataframes with p values
# 
graph_subset_indicator <- F
## Loading the dataframe "df_anpoll_processed"
load(file = paste0(file_path,"anpoll_files/df_anpoll_processed.RData"))
## Loading the dataframe "df_polluter_processed"
load(file = paste0(file_path,"polluter_files/df_polluter_processed.RData"))
load(file = paste0(file_path, "polluter_files/flow_dist_polluter_projected_list.RData"))

## flow_dist_polluter_projected_list <- flow_dist_polluter_projected_cal(file_path)
####################################################
## Loading the analyte-polluter network edgelist "anpoll_edgelist"
load(file = paste0(file_path,"anpoll_files/anpoll_edgelist.RData"))
## Loading the analyte-polluter network edgelist "anpoll_edgelist"
#load(file = paste0(file_path,"anpoll_files/shortest_path_anpoll_edgelist.RData"))
## Loading the "total_edgelist_character_modified"
load(file = paste0(file_path,"common_files_modified/total_edgelist_character_modified.RData"))
## Loading the "flow_dist_from_list"
load(file = paste0(file_path,"inference/flow_dist_from_list.RData"))
## Loading the "flow_dist_to_list"
load(file = paste0(file_path,"inference/flow_dist_to_list.RData"))
## Loading the projected_nodeIDs_list
load(file = paste0(file_path,"polluter_files/projected_nodeIDs_list.RData"))

print("After Loading")


## ## Loading the distances of polluters to projected nodes
## load(file = paste0(file_path,"polluter_files/flow_dist_polluter_projected_list.RData"))

## ## Reading the spill_data_processed to get the affected water body and county
## spill_data_processed <- read.csv(file = paste0(file_path,"Spill_data_processed.csv"), stringsAsFactors = F)
## spill_data_processed$Date <- as.Date(spill_data_processed$Date,format = "%Y-%m-%d")
## str(spill_data_processed)
## str(df_polluter_processed)
## Adding the county information by left joining df_polluter_processed to spill data processed
## df_polluter_processed_appended <- left_join(x = df_polluter_processed, y = spill_data_processed, by = c("lon"="lon", "lat"="lat", "date" = "Date"))
## str(df_polluter_processed_appended)
## save(df_polluter_processed_appended, file = paste0(file_path,"polluter_files/df_polluter_processed_appended.RData"))

## df_polluter_processed_appended[nrow(df_polluter_processed_appended)+1, ] = df_polluter_processed_appended[nrow(df_polluter_processed_appended), ]

df_polluter_processed_appended  <- df_polluter_processed
df_polluter_processed_appended$county  <- 'NA'
####################################################
## Defining the downstream_threshold_dist_km lower and upper
#downstream_threshold_dist_km<-matrix(c(0,10,10,50),2,2,byrow=TRUE)
df_threshold_dist_km <- data.frame("polluter_intersection" = numeric(), "upstream" = numeric(), "downstream_lower"=numeric(), "downstream_upper" = numeric())

df_threshold_dist_km[2,] <- c(5, 5, 0, 10)
df_threshold_dist_km[1,] <- c(45, 5, 10, 50)
print("Before for loop")
## Loop to get test results for spills for downstream close and far
for (j in 1:nrow(df_threshold_dist_km)){
    polluter_test_matrix <- data.frame(matrix(NA, nrow(df_polluter_processed), 15))
    upstream_downstream_obs_list<-list()
    ## Loop to get test results for spills in rank_subgraph
    for(i in 1:6){
        print("Before polluter_test_dist_time()")
        ## for(i in 1:nrow(df_polluter_processed_appended)){
        if(graph_subset_indicator){
            test_result <- polluter_test_dist_time(df_polluter = df_polluter_processed, polluter_lon = df_polluter_processed$lon_mapped[i], polluter_lat = df_polluter_processed$lat_mapped[i], polluter_projected_dist_km = df_threshold_dist_km[j,"polluter_intersection"],upstream_threshold_dist_km = df_threshold_dist_km[j,"upstream"], downstream_threshold_lower_dist_km = df_threshold_dist_km[j,"downstream_lower"], downstream_threshold_upper_dist_km = df_threshold_dist_km[j,"downstream_upper"],date_start=as.Date("1920-01-01"), date_end=as.Date("2017-12-31"), spill_date = df_polluter_processed$date[i], file_path = file_path)
        }else{
            test_result <- polluter_test_dist_time(df_polluter = df_polluter_processed, polluter_lon = df_polluter_processed$lon_mapped[i], polluter_lat = df_polluter_processed$lat_mapped[i], polluter_projected_dist_km = df_threshold_dist_km[j,"polluter_intersection"],upstream_threshold_dist_km = df_threshold_dist_km[j,"upstream"], downstream_threshold_lower_dist_km = df_threshold_dist_km[j,"downstream_lower"], downstream_threshold_upper_dist_km = df_threshold_dist_km[j,"downstream_upper"], date_start=as.Date("1920-01-01"), date_end=as.Date("2017-12-31"), spill_date = df_polluter_processed$date[i], file_path = file_path)
        }
        print("Inside For loop")
        print(i)
                                        #assign("last.warning", NULL, envir = baseenv())
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
    
}
