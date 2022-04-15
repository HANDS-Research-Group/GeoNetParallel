########################################################################################################
file_path<-"/Users/Amal/Box Sync/PSU/Geoscience_Research/River network project/Analysis for Spill Paper/To_Alex/Cl_spill_whole/revision_analysis/3km/"

source(file = paste0(file_path, "code/BaSu_network_v16_modular_functions.R"))

########################################################################################################
## Defining a function for FDR control using BH procedure, Inputs are 1) a vector of p values and 2) significance level alpha (usually set as 0.1). Outputs are 1) a vector of decisions (reject null hypothesis is 1) and 2) a flag indicating if we made atleast one discovery
fdr_decision_cal <- function(p_val,alpha){
  p_val_ordered <- sort(p_val)
  fdr_decision <- rep(0,length(p_val_ordered))
  flag_atleast_1_discovery <- 0
  for (i in length(p_val_ordered):1){
    if(p_val_ordered[i] <= ((i/length(p_val_ordered))*alpha)){
      i_max <- i
      flag_atleast_1_discovery <- 1
      fdr_decision[order(p_val)[1:i_max]] <- 1
      break()
    }
  }
  return(list(fdr_decision, flag_atleast_1_discovery))
}

## Defining a function to genereate final test result dataframe. Inputs are 1) polluter test matrix, 2) significance level alpha (usually set as 0.1) and 3) file_path. Output is the final dataframe with starred p values with or withour fdr control (depending on whether the number of tests is greater than or equal to 15)
fdr_analysis_wrapper <- function(polluter_test_matrix, alpha, file_path){
  
  ## Loading the dataframe df_polluter_processed
  load(file = paste0(file_path,"polluter_files/df_polluter_processed.RData"))
  
  ####################################################
  p_val_rounded_fdr_decision_generator <- function(p_val_mat, alpha, file_path){
    ## FDR Analysis and saving the results
    ## Getting the row IDs in polluter_test_matrix for which we have some test results
    test_pass_ids <- which(p_val_mat[,3] == 1)
    ## Loading the dataframe df_polluter_processed
    load(file = paste0(file_path,"polluter_files/df_polluter_processed.RData"))
    
    p_values_t_test_one_sided_rounded <- round(p_val_mat[,1],3)
    p_values_wilcox_test_one_sided_rounded <- round(p_val_mat[,2],3)
    
    if(length(test_pass_ids) >= 15){
      ## Initializing the fdr decision test vectors as NA
      fdr_decision_t_test_1_sided <- rep(NA_integer_,nrow(p_val_mat))
      fdr_decision_wilcox_test_1_sided <- rep(NA_integer_,nrow(p_val_mat))
      
      ## Getting the FDR decision for Wilcoxon two sided and one side tests
      fdr_decision_t_test_1_sided[test_pass_ids] <- fdr_decision_cal(p_val = p_val_mat[test_pass_ids,1], alpha = alpha)[[1]]
      fdr_decision_wilcox_test_1_sided[test_pass_ids] <- fdr_decision_cal(p_val = p_val_mat[test_pass_ids,2],alpha = alpha)[[1]]
      p_val_mat_output <- cbind("p_values_t_test" = p_values_t_test_one_sided_rounded, "p_values_wilcox_test" =  p_values_wilcox_test_one_sided_rounded, "fdr_decision_t_test" = fdr_decision_t_test_1_sided, "fdr_decision_wilcox_test" = fdr_decision_wilcox_test_1_sided)
    }else{
      p_val_mat_output <- cbind("p_values_t_test" = p_values_t_test_one_sided_rounded, "p_values_wilcox_test" = p_values_wilcox_test_one_sided_rounded, "fdr_decision_t_test" = NA_integer_, "fdr_decision_wilcox_test" = NA_integer_)
    }
    return(p_val_mat_output)
  }
  
  ####################################################
  p_val_up <- p_val_rounded_fdr_decision_generator(p_val_mat = polluter_test_matrix[,c(7,8,9)], alpha = alpha, file_path = file_path)
  p_val_down <- p_val_rounded_fdr_decision_generator(p_val_mat = polluter_test_matrix[,c(10,11,12)], alpha = alpha, file_path = file_path)
  p_val_updown <- p_val_rounded_fdr_decision_generator(p_val_mat = polluter_test_matrix[,c(13,14,15)], alpha = alpha, file_path = file_path)
  
  ## Getting the final indicators for three versions with varying degree of conservativeness
  version_1_pass_mean <- as.numeric((p_val_up[,"p_values_t_test"] > 0.05) & (p_val_down[,"p_values_t_test"] <= 0.05) & (p_val_updown[,"p_values_t_test"] <= 0.05))
  version_1_pass_median <- as.numeric((p_val_up[,"p_values_wilcox_test"] > 0.05) & (p_val_down[,"p_values_wilcox_test"] <= 0.05) & (p_val_updown[,"p_values_wilcox_test"] <= 0.05))
  
  version_2_pass_mean <- as.numeric(p_val_down[,"fdr_decision_t_test"] & p_val_updown[,"fdr_decision_t_test"])
  version_2_pass_median <- as.numeric(p_val_down[,"fdr_decision_wilcox_test"] & p_val_updown[,"fdr_decision_wilcox_test"])
  
  version_3_pass_mean <- as.numeric(p_val_updown[,"fdr_decision_t_test"])
  version_3_pass_median <- as.numeric(p_val_updown[,"fdr_decision_wilcox_test"])
  
  ####################################################
  ## Getting the final dataframe df polluter test mean
  
  df_polluter_test_mean <- data.frame(polluter_test_matrix[,"lon"], polluter_test_matrix[,"lat"], polluter_test_matrix[,"lon_mapped"], polluter_test_matrix[,"lat_mapped"], as.character(df_polluter_processed[,"date"]), round(polluter_test_matrix[,c(1,2)],2),polluter_test_matrix[,c(5,6)], p_val_up[,"p_values_t_test"], p_val_down[,"p_values_t_test"], p_val_updown[,"p_values_t_test"], version_1_pass_mean, version_2_pass_mean, version_3_pass_mean, stringsAsFactors = F)
  names(df_polluter_test_mean) <- c("Longitude_Original", "Latitude_Original", "Longitude_Mapped", "Latitude_Mapped", "Date", "Upstream Mean (ppb)","Downstream Mean (ppb)","No. of Observations Upstream","No. of Observations Downstream","t test (up) p values","t test (down) p values", "t test (updown) p values", "version_1_decision", "version_2_decision", "version_3_decision")
  
  ## Filtering for volume >= 400 gallons
  #df_polluter_test_mean <- df_polluter_test_mean %>% filter((!is.na(Volume)) & (Volume >= 400))
  
  ####################################################
  ## Getting the final dataframe df polluter test median
  df_polluter_test_median <- data.frame(polluter_test_matrix[,"lon"], polluter_test_matrix[,"lat"], polluter_test_matrix[,"lon_mapped"], polluter_test_matrix[,"lat_mapped"], as.character(df_polluter_processed[,"date"]), round(polluter_test_matrix[,c(3,4)],2),polluter_test_matrix[, c(5,6)], p_val_up[,"p_values_wilcox_test"], p_val_down[,"p_values_wilcox_test"], p_val_updown[,"p_values_wilcox_test"],  version_1_pass_median, version_2_pass_median, version_3_pass_median, stringsAsFactors = F)
  names(df_polluter_test_median) <- c("Longitude_Original", "Latitude_Original", "Longitude_Mapped", "Latitude_Mapped", "Date","Upstream Median (ppb)","Downstream Median (ppb)","No. of Observations Upstream","No. of Observations Downstream","Wilcox (up) p values","Wilcox (down) p values","Wilcox (updown) p values", "version_1_decision", "version_2_decision", "version_3_decision")
  
  ## Filtering for volume >= 400 gallons
  #df_polluter_test_median <- df_polluter_test_median %>% filter((!is.na(Volume)) & (Volume >= 400))
  
  return(list(df_polluter_test_mean, df_polluter_test_median))
}

########################################################################################################
########################################################################################################
## Loading the test results and processing them in a single matrix and final dataframes

## Defining the downstream_threshold_dist_km lower and upper
df_threshold_dist_km<-data.frame("polluter_intersection"=numeric(),"upstream"=numeric(),"downstream_lower"=numeric(),"downstream_upper"=numeric())

df_threshold_dist_km[1,] <- c(5,5,0,10)
df_threshold_dist_km[2,] <- c(45,5,10,50)

## setting the j corresponding to second set of distance parameters in df_threshold_dist_km
j <- 2

## Loading the appended polluter processed with county info.
load(file = paste0(file_path,"polluter_files/df_polluter_processed.RData"))

## Initializing the matrix and observation storage lists
polluter_test_matrix <- data.frame(matrix(NA,nrow(df_polluter_processed),19))
upstream_downstream_obs_list<-list()

#missing_ids<-c()
for (i in 1:nrow(polluter_test_matrix)){
  if(class(try(load(file = paste0(file_path,"inference/test_results_j",j,"/test_result_",i,".RData")),silent = T))=="try-error"){
    missing_ids<-c(missing_ids,i)
  }else{
    print(i)
    load(file = paste0(file_path,"inference/test_results_j",j,"/test_result_",i,".RData"))
    ## Getting mean, median and no. of observations for upstream and downstream in the first 6 columns
    polluter_test_matrix[i,c(1:6)] <- test_result[[3]][c(1:6)]
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
    polluter_test_matrix[i,15] <- as.numeric(test_result[[3]][9])
    
    ## Adding original latitude and longitude
    polluter_test_matrix[i,16] <- df_polluter_processed$lon[i]
    polluter_test_matrix[i,17] <- df_polluter_processed$lat[i]
    
    ## Adding mapped latitude and longitude
    polluter_test_matrix[i,18] <- df_polluter_processed$lon_mapped[i]
    polluter_test_matrix[i,19] <- df_polluter_processed$lat_mapped[i]
    
    ## Initializing and storing the observations for upstream and downstream
    upstream_downstream_obs_list[[i]] <- list()
    upstream_downstream_obs_list[[i]][[1]] <- test_result[[4]]
    upstream_downstream_obs_list[[i]][[2]] <- test_result[[5]]
    #print(i)
    print(polluter_test_matrix[i,])
  }
}
#missing_ids
#save(missing_ids,file=paste0(file_path,"inference/missing_ids_v4.RData"))

colnames(polluter_test_matrix) <- c("upstream mean","downstream mean","upstream median","downstream median","number of obs upstream","number of obs downstream", "t test (up) p value", "Wilcoxon test (up) p value","test result indicator (up)","t test (down) p value", "Wilcoxon test (down) p value","test result indicator (down)","t test (updown) p value", "Wilcoxon test (updown) p value","test result indicator (updown)", "lon", "lat", "lon_mapped", "lat_mapped")

save(polluter_test_matrix,file = paste0(file_path,"inference/polluter_test_matrix_",df_threshold_dist_km[j,"downstream_upper"],".RData"))

#####################################################
## Doing the fdr analysis for >=15 pvalues in each test
load(file = paste0(file_path,"inference/polluter_test_matrix_",df_threshold_dist_km[j,"downstream_upper"],".RData"))

#debug(fdr_analysis_wrapper)
df_polluter_test <- fdr_analysis_wrapper(polluter_test_matrix = polluter_test_matrix, alpha = 0.1, file_path = file_path)

df_polluter_test_mean <- df_polluter_test[[1]]
df_polluter_test_median <- df_polluter_test[[2]]
str(df_polluter_test_mean)
str(df_polluter_test_median)

save(df_polluter_test_mean,file = paste0(file_path,"inference/df_polluter_test_mean_", df_threshold_dist_km[j,"downstream_upper"], ".RData"))
save(df_polluter_test_median,file = paste0(file_path,"inference/df_polluter_test_median_", df_threshold_dist_km[j,"downstream_upper"], ".RData"))

write.csv(df_polluter_test_mean,file = paste0(file_path,"inference/df_polluter_test_mean_", df_threshold_dist_km[j,"downstream_upper"], ".csv"))
write.csv(df_polluter_test_median,file = paste0(file_path,"inference/df_polluter_test_median_", df_threshold_dist_km[j,"downstream_upper"], ".csv"))

str(df_polluter_test_median)

df_polluter_test_median$distance_Tioga_1 <- NA

for(i in 1:nrow(df_polluter_test_median)){
  df_polluter_test_median$distance_Tioga_1[i] <- distm(x = as.matrix(df_polluter_test_median[i,c("Longitude_Original","Latitude_Original")]), y = t(as.matrix(c(-76.965211,41.698031))))
}

str(df_polluter_test_median)

df_polluter_test_median$distance_order <- order(df_polluter_test_median$distance_Tioga_1)

library(dplyr)
# sort the dataframe in R
df_polluter_test_median <- df_polluter_test_median %>% arrange(distance_Tioga_1)

write.csv(df_polluter_test_median,file = paste0(file_path,"inference/df_polluter_test_median_", df_threshold_dist_km[j,"downstream_upper"], ".csv"))


########################################################################################################
# Analysis for validation with Anna's version

# str(df_polluter_test_median)
# 
# ## Left joining with Anna's spills
# df_Anna_spills <- read.csv(file = "/Users/Amal/Box Sync/PSU/Summer 2019/Geoscience_Research/River network project/Analysis for Spill Paper/To_Alex/Anna_spills/v2/List_of_Spills_Sent_to_Amal.csv")
# #df_Anna_spills <- read.csv(file = "/Users/Amal/Box Sync/PSU/Summer 2019/Geoscience_Research/River network project/Analysis for Spill Paper/To_Alex/List_of_Spills_Sent_to_Amal.csv")
# 
# df_Anna_spills$Date_start <- as.Date(df_Anna_spills$Date_start, format = "%m/%d/%y")
# df_Anna_spills$Date_end <- as.Date(df_Anna_spills$Date_end, format = "%m/%d/%y")
# 
# str(df_Anna_spills)
# head(df_Anna_spills)
# 
# df_polluter_test_median$Date <- as.Date(df_polluter_test_median$Date)
# str(df_polluter_test_median)
# 
# #between(as.Date("2010-05-24"), df_Anna_spills$Date_start[1], df_Anna_spills$Date_end[1])
# 
# #df_polluter_test_median_Anna_appended <- left_join(x = df_polluter_test_median, y = df_Anna_spills, by = c("Longitude_Original"="Long", "Latitude_Original" = "Lat", "Date" = "Date.of.Incident"), all.x = TRUE)
# df_polluter_test_median_Anna_appended <- left_join(x = df_polluter_test_median, y = df_Anna_spills, by = c("Longitude_Original"="Long", "Latitude_Original" = "Lat"), all.x = TRUE)
# 
# str(df_polluter_test_median_Anna_appended)
# head(df_polluter_test_median_Anna_appended)
# sum(df_polluter_test_median_Anna_appended$Anna_Spills, na.rm = T)
# for(i in 1:nrow(df_polluter_test_median_Anna_appended)){
#   if(!is.na(df_polluter_test_median_Anna_appended$Anna_Spills[i])){
#     if(df_polluter_test_median_Anna_appended$Anna_Spills[i]==1){
#       if(!between(df_polluter_test_median_Anna_appended$Date[i], df_polluter_test_median_Anna_appended$Date_start[i], df_polluter_test_median_Anna_appended$Date_end[i])){
#         df_polluter_test_median_Anna_appended$Anna_Spills[i] <- NA
#         df_polluter_test_median_Anna_appended$Anna_Significant_Spills[i] <- NA
#       }
#     }
#   }
# }
# sum(df_polluter_test_median_Anna_appended$Anna_Spills, na.rm = T)
# 
# write.csv(df_polluter_test_median_Anna_appended,file = paste0(file_path,"inference/df_polluter_test_median_", df_threshold_dist_km[j,"downstream_upper"], "_Anna_appended.csv"))

########################################################################################################
# for (i in 1:length(missing_ids)){
#   test_result <- polluter_test_dist_time(df_polluter = df_polluter_processed, polluter_lon = df_polluter_processed$lon_mapped[missing_ids[i]], polluter_lat = df_polluter_processed$lat_mapped[missing_ids[i]], polluter_projected_dist_km = df_threshold_dist_km[j,"polluter_intersection"],upstream_threshold_dist_km = df_threshold_dist_km[j,"upstream"], downstream_threshold_lower_dist_km = df_threshold_dist_km[j,"downstream_lower"], downstream_threshold_upper_dist_km = df_threshold_dist_km[j,"downstream_upper"], date_start=as.Date("1920-01-01"), date_end=as.Date("2017-12-31"), spill_date = df_polluter_processed$date[missing_ids[i]], file_path = file_path)
#   print(i)
# }
# 
# i=1
# debug(polluter_test_dist_time)
# test_result <- polluter_test_dist_time(df_polluter = df_polluter_processed, polluter_lon = df_polluter_processed$lon_mapped[missing_ids[i]], polluter_lat = df_polluter_processed$lat_mapped[missing_ids[i]], polluter_projected_dist_km = df_threshold_dist_km[j,"polluter_intersection"],upstream_threshold_dist_km = df_threshold_dist_km[j,"upstream"], downstream_threshold_lower_dist_km = df_threshold_dist_km[j,"downstream_lower"], downstream_threshold_upper_dist_km = df_threshold_dist_km[j,"downstream_upper"], date_start=as.Date("1920-01-01"), date_end=as.Date("2017-12-31"), spill_date = df_polluter_processed$date[missing_ids[i]], file_path = file_path)

########################################################################################################

# df_spills <- read.csv(file = "/Users/Amal/Box Sync/PSU/Summer 2019/Geoscience_Research/River network project/Analysis for Spill Paper/To_Alex/spill_data_processed_with_volume/Spill_data_processed_v_628.csv")
# 
# str(df_spills)
# head(df_spills)
# nrow(df_spills)
# 
# df_spills %>% select(Volume_Spill) %>% filter(!is.na(Volume_Spill)) %>% filter(Volume_Spill>=400)
# 
# load(file = paste0(file_path, "inference/df_polluter_test_mean_50.RData"))
# str(df_polluter_test_mean)

load(file = paste0(file_path,"inference/df_polluter_test_median_", df_threshold_dist_km[j,"downstream_upper"], ".RData"))

str(df_polluter_test_median)
df_polluter_test_median[16,]

distm(x = as.matrix(df_polluter_test_median[16,c("Longitude_Original","Latitude_Original")]), y = as.matrix(df_polluter_test_median[16,c("Longitude_Mapped","Latitude_Mapped")]))


