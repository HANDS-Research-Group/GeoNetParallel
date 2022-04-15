library(tidyverse)

## Update the path to the folder
file_path <- "/home/rrpatil/GeoNet_2021_SC/"

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

flow_dist_polluter_projected_list <- flow_dist_polluter_projected_cal(file_path)
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
## ## Loading the distances of polluters to projected nodes
## load(file = paste0(file_path,"polluter_files/flow_dist_polluter_projected_list.RData"))

## ## Reading the spill_data_processed to get the affected water body and county
## spill_data_processed <- read.csv(file = paste0(file_path,"Spill_data_processed.csv"), stringsAsFactors = F)
## spill_data_processed$Date <- as.Date(spill_data_processed$Date,format = "%Y-%m-%d")
## str(spill_data_processed)
## str(df_polluter_processed)
## ## Adding the county information by left joining df_polluter_processed to spill data processed
## df_polluter_processed_appended <- left_join(x = df_polluter_processed, y = spill_data_processed, by = c("lon"="lon", "lat"="lat", "date" = "Date"))
## str(df_polluter_processed_appended)
## save(df_polluter_processed_appended, file = paste0(file_path,"polluter_files/df_polluter_processed_appended.RData"))

df_polluter_processed_appended  <- df_polluter_processed

df_threshold_dist_km <- data.frame("polluter_intersection" = numeric(), "upstream" = numeric(), "downstream_lower"=numeric(), "downstream_upper" = numeric())

## df_threshold_dist_km[1,] <- c(5, 5, 0, 10)
df_threshold_dist_km[1,] <- c(45, 5, 0, 50)


# The following function peforms the inference on all the polluter sites already loaded with the GeoNet.
all_polluter_inference_test <- function(df_threshold_km) {
    ####################################################
    ## Defining the downstream_threshold_dist_km lower and upper
                                        #downstream_threshold_dist_km<-matrix(c(0,10,10,50),2,2,byrow=TRUE)

    ## Loop to get test results for spills for downstream close and far
    for (j in 1:nrow(df_threshold_dist_km)){
        polluter_test_matrix <- data.frame(matrix(NA, nrow(df_polluter_processed), 15))
        upstream_downstream_obs_list<-list()
        ## Loop to get test results for spills in rank_subgraph
        for(i in 1:nrow(df_polluter_processed_appended)){
            if(graph_subset_indicator){
                test_result <- polluter_test_dist_time(df_polluter = df_polluter_processed, polluter_lon = df_polluter_processed$lon_mapped[i], polluter_lat = df_polluter_processed$lat_mapped[i], polluter_projected_dist_km = df_threshold_dist_km[j,"polluter_intersection"],upstream_threshold_dist_km = df_threshold_dist_km[j,"upstream"], downstream_threshold_lower_dist_km = df_threshold_dist_km[j,"downstream_lower"], downstream_threshold_upper_dist_km = df_threshold_dist_km[j,"downstream_upper"],date_start=as.Date("1920-01-01"), date_end=as.Date("2017-12-31"), spill_date = df_polluter_processed$date[i], file_path = file_path)
            }else{
                test_result <- polluter_test_dist_time(df_polluter = df_polluter_processed, polluter_lon = df_polluter_processed$lon_mapped[i], polluter_lat = df_polluter_processed$lat_mapped[i], polluter_projected_dist_km = df_threshold_dist_km[j,"polluter_intersection"],upstream_threshold_dist_km = df_threshold_dist_km[j,"upstream"], downstream_threshold_lower_dist_km = df_threshold_dist_km[j,"downstream_lower"], downstream_threshold_upper_dist_km = df_threshold_dist_km[j,"downstream_upper"], date_start=as.Date("1920-01-01"), date_end=as.Date("2017-12-31"), spill_date = df_polluter_processed$date[i], file_path = file_path)
            }
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
}

# The following function takes in as input the inference output on all
# the polluter sites already present in GeoNet and a dataframe with
# new temporal polluter information and the threshold distance in km
# to produce new inference result
add_temporal_polluter_inference_test <- function(df_polluter_processed, df_threshold_km) {

    df_polluter_processed_appended  <- df_polluter_processed
    for (j in 1:nrow(df_threshold_dist_km)){
        polluter_test_matrix_to_append <- data.frame(matrix(NA, nrow(df_polluter_processed), 15))
        upstream_downstream_obs_list<-list()
        ## Loop to get test results for spills in rank_subgraph
        for(i in 1:nrow(df_polluter_processed_appended)){
            if(graph_subset_indicator){
                test_result <- polluter_test_dist_time(df_polluter = df_polluter_processed, polluter_lon = df_polluter_processed$lon_mapped[i], polluter_lat = df_polluter_processed$lat_mapped[i], polluter_projected_dist_km = df_threshold_dist_km[j,"polluter_intersection"],upstream_threshold_dist_km = df_threshold_dist_km[j,"upstream"], downstream_threshold_lower_dist_km = df_threshold_dist_km[j,"downstream_lower"], downstream_threshold_upper_dist_km = df_threshold_dist_km[j,"downstream_upper"],date_start=as.Date("1920-01-01"), date_end=as.Date("2017-12-31"), spill_date = df_polluter_processed$date[i], file_path = file_path)
            }else{
                test_result <- polluter_test_dist_time(df_polluter = df_polluter_processed, polluter_lon = df_polluter_processed$lon_mapped[i], polluter_lat = df_polluter_processed$lat_mapped[i], polluter_projected_dist_km = df_threshold_dist_km[j,"polluter_intersection"],upstream_threshold_dist_km = df_threshold_dist_km[j,"upstream"], downstream_threshold_lower_dist_km = df_threshold_dist_km[j,"downstream_lower"], downstream_threshold_upper_dist_km = df_threshold_dist_km[j,"downstream_upper"], date_start=as.Date("1920-01-01"), date_end=as.Date("2017-12-31"), spill_date = df_polluter_processed$date[i], file_path = file_path)
            }
            print(i)
                                        #assign("last.warning", NULL, envir = baseenv())
            ## Getting mean, median and no. of observations for upstream and downstream in the first 6 columns
            polluter_test_matrix_to_append[i, c(1:6)] <- test_result[[3]][c(1:6)]
            ## Storing the p values and test indicator for upstream temporal test
            if(as.numeric(test_result[[1]][3]) == 1){
                polluter_test_matrix_to_append[i,7] <- as.numeric(test_result[[1]][1])
                polluter_test_matrix_to_append[i,8] <- as.numeric(test_result[[1]][2])
            }
            polluter_test_matrix_to_append[i,9] <- as.numeric(test_result[[1]][3])
            ## Storing the p values and test indicator for downstream temporal test
            if(as.numeric(test_result[[2]][3]) == 1){
                polluter_test_matrix_to_append[i,10] <- as.numeric(test_result[[2]][1])
                polluter_test_matrix_to_append[i,11] <- as.numeric(test_result[[2]][2])
            }
            polluter_test_matrix_to_append[i,12] <- as.numeric(test_result[[2]][3])
            ## Storing the p values and test indicator for final upstream vs. downstream spatio-temporal test
            if(as.numeric(test_result[[3]][9]) == 1){
                polluter_test_matrix_to_append[i,13] <- as.numeric(test_result[[3]][7])
                polluter_test_matrix_to_append[i,14] <- as.numeric(test_result[[3]][8])
            }
            polluter_test_matrix_to_append[i,15]<-as.numeric(test_result[[3]][9])
            
                                        # Lat/Lon information
            polluter_test_matrix_to_append[i,16]<-df_polluter_processed_appended$lon[i]
            polluter_test_matrix_to_append[i,17]<-df_polluter_processed_appended$lat[i]
            polluter_test_matrix_to_append[i,18]<-df_polluter_processed_appended$lon_mapped[i]
            polluter_test_matrix_to_append[i,19]<-df_polluter_processed_appended$lat_mapped[i]
            
            ## Initializing and storing the observations for upstream and downstream
            upstream_downstream_obs_list[[i]] <- list()
            upstream_downstream_obs_list[[i]][[1]] <- test_result[[4]]
            upstream_downstream_obs_list[[i]][[2]] <- test_result[[5]]
                                        #print(i)
            print(polluter_test_matrix_to_append[i,])
        }
        
        colnames(polluter_test_matrix_to_append) <- c("upstream mean","downstream mean","upstream median","downstream median","number of obs upstream","number of obs downstream", "t test (up) p value", "Wilcoxon test (up) p value","test result indicator (up)","t test (down) p value", "Wilcoxon test (down) p value","test result indicator (down)","t test (updown) p value", "Wilcoxon test (updown) p value","test result indicator (updown)", "lon", "lat", "lon_mapped", "lat_mapped")
        load(file = paste0(file_path,"inference/polluter_test_matrix_",df_threshold_dist_km[j,"downstream_upper"],".RData"))
        polluter_test_matrix = rbind(polluter_test_matrix, polluter_test_matrix_to_append)
        save(polluter_test_matrix, file = paste0(file_path, "inference/polluter_test_matrix_",df_threshold_dist_km[j,"downstream_upper"],".RData"))
    }
    
}


## Function takes the dataframe with new temporal analyte information
## and re-runs the inference for all polluter sites
add_temporal_analyte_inference_test <- function(df_analyte_to_append, df_threshold_km) {
    load(file = paste0(file_path, "analyte_files/df_anpoll_processed.RData"))    
    df_anpoll_processed  <- rbind(df_anpoll_processed, df_analyte_to_append)
    df_list_analyte_obj <- df_list_analyte_generator(df_analyte_processed = df_anpoll_processed[which(df_anpoll_processed$anpoll_indicator=="analyte"),],output_path = file_path)
    all_polluter_inference_test(df_threshold_km)
}

## The following function to add both analyte and polluter together at a time.
add_temporal_analyte_polluter_inference_test <- function(df_analyte_to_append, df_polluter_to_append, df_threshold_km) {
    add_temporal_polluter_inference_test(df_polluter_to_append, df_threshold_km)
    add_temporal_analyte_inference_test(df_analyte_to_append, df_threshold_km)
}



### Test Code
## Run test for all polluters
## all_polluter_inference_test()

## Add new polluter existing location and run test
load(file = paste0(file_path,"polluter_files/df_polluter_to_append.RData"))

df_polluter_to_append <- df_polluter_processed[which(df_polluter_processed$lon==df_polluter_preprocessed$lon[1]),]
df_polluter_to_append <- max(df_polluter_processed$date + 1)
# df_polluter_to_append has id, date, lat, lon 
add_temporal_polluter_inference_test(df_polluter_to_append, df_threshold_km)


## ## ## Add new analyte at existing location and run test
## load(file = paste0(file_path, "analyte_files/df_anpoll_processed.RData"))    
## df_analyte_to_append <- df_anpoll_processed[which(df_anpoll_processed$lon==df_anpoll_processed$lon[1] &&
##                                                         df_anpoll_processed$lat==df_anpoll_processed$lat[1]), ]

## df_analyte_to_append$date <- max(df_anpoll_processed$date+1)
## # df_analyte_to_append has id, date, lat, lon, conc
## add_temporal_analyte_inference_test(df_analyte_to_append, df_threshold_km)

