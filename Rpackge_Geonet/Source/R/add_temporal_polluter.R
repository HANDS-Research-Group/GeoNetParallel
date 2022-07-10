#' Some title
#' @description A description of this function with inline latex \eqn{2^3=8}
#' test for displayed equations \deqn{\frac{1}{2}\sum_{i=1}^n\log(x+y)}.
#' We can also use links, tables, lists and Character formattings such as bold and italic.
#' @param path_total description of this parameter
#' @param df_polluter_to_append_file_path description of this parameter
#' @param upstream_thresh description of this parameter
#' @param downstream_lower_thresh description of this parameter
#' @param downstream_upper_thresh description of this parameter
#' @param num_cores description of this parameter
#' @param permission description of this parameter
#'
#' @return description of the return value
#' @export
#'
#' @useDynLib Rpackage0621
#' @importFrom Rcpp sourceCpp
#' @importFrom foreach %do%
#' @importFrom foreach %dopar%

add_temporal_polluter  <- function(path_total, df_polluter_to_append_file_path,
                                   upstream_thresh = 5,
                                   downstream_lower_thresh = 0,
                                   downstream_upper_thresh = 50,
                                   num_cores = 1, permission = "no") {
  if (!(permission == "yes")){
    print("You need to explicit set this variable to make sure that you know what will happen.")
    return (0)
  }
  file_path = path_total
  n_chunks = num_cores
  ##################################################
  df_threshold_dist_km <- data.frame("polluter_intersection" = numeric(),
                                     "upstream" = numeric(), "downstream_lower" = numeric(),
                                     "downstream_upper" = numeric())
  df_threshold_dist_km <- data.frame("polluter_intersection" = numeric(),"upstream" = numeric(), "downstream_lower" = numeric(),"downstream_upper" = numeric())

  df_threshold_dist_km[1,] <- c(downstream_upper_thresh-5, upstream_thresh, downstream_lower_thresh, downstream_upper_thresh)

  df_polluter_raw<-utils::read.csv(file = df_polluter_to_append_file_path)

  ## Preprocessing the analyte dataframe
  df_polluter_preprocessed<-df_polluter_raw[,c(2,3)]
  df_polluter_preprocessed$date<-as.Date(df_polluter_raw$date,format="%m/%d/%Y")

  load(file=paste0(file_path,"polluter_files/df_polluter_processed.RData"))

  df_polluter_preprocessed$lon_mapped <- NA
  df_polluter_preprocessed$lat_mapped <- NA
  df_polluter_preprocessed$nodeID <- NA
  df_polluter_preprocessed$conc <- NA
  df_polluter_preprocessed$anpoll_indicator <- "polluter"
  df_polluter_preprocessed$lon_mapped <- df_polluter_processed$lon_mapped[which(df_polluter_preprocessed$lon==df_polluter_processed$lon &
                                                                                  df_polluter_preprocessed$lat==df_polluter_processed$lat)]
  df_polluter_preprocessed$lat_mapped <- df_polluter_processed$lat_mapped[which(df_polluter_preprocessed$lat==df_polluter_processed$lat &
                                                                                  df_polluter_preprocessed$lon==df_polluter_processed$lon)]
  df_polluter_preprocessed$nodeID <- df_polluter_processed$nodeID[which(df_polluter_preprocessed$lon==df_polluter_processed$lon &
                                                                          df_polluter_preprocessed$lat==df_polluter_processed$lat)]

  graph_subset_indicator <- F
  df_polluter_processed_appended  <- df_polluter_preprocessed
  for (j in 1:nrow(df_threshold_dist_km)){
    polluter_test_matrix_to_append <- data.frame(matrix(NA, nrow(df_polluter_preprocessed), 19))
    upstream_downstream_obs_list<-list()
    ## Loop to get test results for spills in rank_subgraph
    for(i in 1:nrow(df_polluter_preprocessed)){
      if(graph_subset_indicator){
        test_result <- polluter_test_dist_time(df_polluter = df_polluter_preprocessed, polluter_lon = df_polluter_preprocessed$lon_mapped[i], polluter_lat = df_polluter_preprocessed$lat_mapped[i], polluter_projected_dist_km = df_threshold_dist_km[j,"polluter_intersection"],upstream_threshold_dist_km = df_threshold_dist_km[j,"upstream"], downstream_threshold_lower_dist_km = df_threshold_dist_km[j,"downstream_lower"], downstream_threshold_upper_dist_km = df_threshold_dist_km[j,"downstream_upper"],date_start=as.Date("1920-01-01"), date_end=as.Date("2017-12-31"), spill_date = df_polluter_processed$date[i], file_path = file_path)
      }else{
        test_result <- polluter_test_dist_time(df_polluter = df_polluter_preprocessed, polluter_lon = df_polluter_preprocessed$lon_mapped[i], polluter_lat = df_polluter_preprocessed$lat_mapped[i], polluter_projected_dist_km = df_threshold_dist_km[j,"polluter_intersection"],upstream_threshold_dist_km = df_threshold_dist_km[j,"upstream"], downstream_threshold_lower_dist_km = df_threshold_dist_km[j,"downstream_lower"], downstream_threshold_upper_dist_km = df_threshold_dist_km[j,"downstream_upper"], date_start=as.Date("1920-01-01"), date_end=as.Date("2017-12-31"), spill_date = df_polluter_processed$date[i], file_path = file_path)
      }
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
    }

    colnames(polluter_test_matrix_to_append) <- c("upstream mean","downstream mean","upstream median","downstream median","number of obs upstream","number of obs downstream", "t test (up) p value", "Wilcoxon test (up) p value","test result indicator (up)","t test (down) p value", "Wilcoxon test (down) p value","test result indicator (down)","t test (updown) p value", "Wilcoxon test (updown) p value","test result indicator (updown)", "lon", "lat", "lon_mapped", "lat_mapped")
    load(file = paste0(file_path,"inference/polluter_test_matrix_",df_threshold_dist_km[j,"downstream_upper"],".RData"))
    polluter_test_matrix = rbind(polluter_test_matrix, polluter_test_matrix_to_append)
    save(polluter_test_matrix, file = paste0(file_path, "inference/polluter_test_matrix_",df_threshold_dist_km[j,"downstream_upper"],".RData"))
  }

  ####################
  ## FDR
  ####################

  df_polluter_processed_appended <- rbind(df_polluter_processed, df_polluter_processed_appended)
  df_polluter_processed <- df_polluter_processed_appended
  save(df_polluter_processed, file = paste0(file_path, "polluter_files/df_polluter_processed.RData"))

  df_polluter_processed_appended$County  <-  NA
  ## df_polluter_processed_appended[nrow(df_polluter_processed_appended)+1,]  <- df_polluter_processed_appended[nrow(df_polluter_processed_appended),]
  for(j in 1:nrow(df_threshold_dist_km)) {
    ## for (j in 2:nrow(df_threshold_dist_km)){

    #####################################################
    ## Doing the fdr analysis for >=15 pvalues in each test
    load(file = paste0(file_path,"inference/polluter_test_matrix_",df_threshold_dist_km[j,"downstream_upper"],".RData"))
    df_polluter_test <- fdr_analysis_wrapper(polluter_test_matrix = polluter_test_matrix,alpha = 0.1,county=df_polluter_processed_appended$County, file_path = file_path)
    df_polluter_test_mean <- df_polluter_test[[1]]
    df_polluter_test_median <- df_polluter_test[[2]]
    save(df_polluter_test_mean, file = paste0(file_path,"inference/df_polluter_test_mean_", df_threshold_dist_km[j, "downstream_upper"], ".RData"))
    save(df_polluter_test_median, file = paste0(file_path,"inference/df_polluter_test_median_", df_threshold_dist_km[j, "downstream_upper"], ".RData"))
  }

}
