#' Some title
#' @description A description of this function with inline latex \eqn{2^3=8}
#' test for displayed equations \deqn{\frac{1}{2}\sum_{i=1}^n\log(x+y)}.
#' We can also use links, tables, lists and Character formattings such as bold and italic.
#' @param path_total description of this parameter
#' @param df_analyte_to_append_filepath description of this parameter
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

add_temporal_analyte_polluter <- function(path_total, df_analyte_to_append_filepath, df_polluter_to_append_file_path,
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
  add_temporal_polluter(path_total, df_polluter_to_append_file_path,
                        upstream_thresh,
                        downstream_lower_thresh,
                        downstream_upper_thresh,
                        num_cores, permission)
  add_temporal_analyte(path_total, df_analyte_to_append_filepath,
                       upstream_thresh,
                       downstream_lower_thresh,
                       downstream_upper_thresh,
                       num_cores, permission)
}
