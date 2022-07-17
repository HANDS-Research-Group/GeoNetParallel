#' Add new water chemistry data and polluting events at new locations
#' @description This function can be used to add new water chemistry data and polluting events at new locations and run the statistical inference for all the polluting events. This function can only be used after “whole_pipe” is executed. Execution of this function will overwrite the previously generated output. Therefore, it is advised to preserve a copy of the old outputs in case you need it.
#' @param path_total This is the same root file path as that of the function “whole_pipe”. The previous output will be overwritten so backup is recommended.
#' @param df_analyte_to_append_filepath This is the full path to the csv file containing new water chemistry data for the given stream. The file formatting requirement is the same as that of “filename_chem” in “whole_pipe”.
#' @param df_polluter_to_append_file_path This is the full path to the csv file containing new pollution site data for the given stream. The file formatting requirement is the same as that of “filename_poll” in “whole_pipe”.
#' @param upstream_thresh The same parameter as that in "whole_pipe".
#' @param downstream_lower_thresh The same parameter as that in "whole_pipe".
#' @param downstream_upper_thresh The same parameter as that in "whole_pipe".
#' @param num_cores The same parameter as that in "whole_pipe".
#' @param permission The same parameter as that in "whole_pipe".
#' @return The same return files as those in "whole_pipe".
#' @note This function assumes “whole_pipe” function has already been run and the outputs generated from it lie in “path_total” folder. Please refer to our \href{https://github.com/HANDS-Research-Group/StreamNet/tree/main/Rpackage_StreamNet}{website} for examples.
#' @export
#'
#' @useDynLib StreamNet
#' @importFrom Rcpp sourceCpp
#' @importFrom foreach %do%
#' @importFrom foreach %dopar%

add_new_analyte_polluter <- function(path_total, df_analyte_to_append_filepath, df_polluter_to_append_file_path,
                                     upstream_thresh = 5,
                                     downstream_lower_thresh = 0,
                                     downstream_upper_thresh = 50,
                                     num_cores = 1, permission = "no") {
  if (!(permission == "yes")){
    print("You need to explicit set the variable permission to make sure that you know what will happen. Please refer to the help page for function whole_pipe.")
    return (0)
  }
  file_path = path_total
  n_chunks = num_cores
  ##################################################
  add_new_polluter(path_total, df_polluter_to_append_file_path,
                   upstream_thresh,
                   downstream_lower_thresh,
                   downstream_upper_thresh,
                   num_cores, permission)
  add_new_analyte(path_total, df_analyte_to_append_filepath,
                  upstream_thresh,
                  downstream_lower_thresh,
                  downstream_upper_thresh,
                  num_cores, permission)
}
