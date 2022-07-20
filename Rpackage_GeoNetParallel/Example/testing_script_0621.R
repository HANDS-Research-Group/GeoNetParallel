setwd("/home/handslab/Desktop/zhong_zheng/test_envir_medium_all_0621")
library(GeoNetParallel)
############test whole pipeline################################
path_total = "/home/handslab/Desktop/zhong_zheng/test_envir_medium_all_0621/test_envir_whole_pipe/"
filename_chem = "Water_Chemistry_R_Package.csv"
filename_poll = "Pollution_Site_R_Package.csv"
num_cores = 26
permission = "yes"
try(whole_pipe(path_total, filename_chem, filename_poll, num_cores, permission))
print("task 1 finished")
###########test add new analy#################
path_total = "/home/handslab/Desktop/zhong_zheng/test_envir_medium_all_0621/test_envir_add_new_anl/"
filename_chem = "Water_Chemistry_R_Package_missing.csv"
filename_poll = "Pollution_Site_R_Package.csv"
df_analyte_to_append_filepath = "/home/handslab/Desktop/zhong_zheng/test_envir_medium_all_0621/test_envir_add_new_anl/data/Water_chemistry_new.csv"
num_cores = 26
permission = "yes"
try(whole_pipe(path_total, filename_chem, filename_poll, num_cores, permission))
try(add_new_analyte(path_total, df_analyte_to_append_filepath, num_cores, permission))
print("task 2 finished")
########test add temp anl######################
path_total = "/home/handslab/Desktop/zhong_zheng/test_envir_medium_all_0621/test_envir_add_temp_anl/"
filename_chem = "Water_Chemistry_R_Package_missing_temporal.csv"
filename_poll = "Pollution_Site_R_Package.csv"
df_analyte_to_append_filepath = "/home/handslab/Desktop/zhong_zheng/test_envir_medium_all_0621/test_envir_add_temp_anl/data/Water_chemistry_temporal.csv"
num_cores = 26
permission = "yes"
try(whole_pipe(path_total, filename_chem, filename_poll, num_cores, permission))
try(add_temporal_analyte(path_total, df_analyte_to_append_filepath, num_cores, permission))
print("task 3 finished")
