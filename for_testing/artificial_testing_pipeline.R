# Title     : TODO
# Objective : TODO
# Created by: fabio
# Created on: 02/03/2021

# rm(list = ls())

library(pcalg)
library(infotheo)
library(BDgraph)
library(sbgcop)

library(ConfigParser)
# from demo_copulaFactorPC.R

source('./for_testing/artificial_generator_data_and_factor.R')

source( './R_code/backend_galaina_final_part.R')

# if (sys.nframe() == 0) {
# main <- function()
# {
tmp_path_to_basic_config <- '/Users/fabio/projects/aggressotype/code/galaina/user_data/test/AsiaSmall/config_5/config.ini'
# TODO - make a config.ini paying attention at:
# input_path_directed_edges_blacklist

dataset_n <-3
fileConn <- NULL

config <- ConfigParser$new(Sys.getenv(), optionxform = identity)
config$read(filepath = tmp_path_to_basic_config)
load_copula_factor_functions(config)

multi_out_ls <- multi_artificial_generator_data_and_factor(dataset_n = dataset_n )

i_run <- 1
tmp_Y <- multi_out_ls[[i_run]][['for_inference']][['Y']]
tmp_Lambda <- multi_out_ls[[i_run]][['for_inference']][['Lambda']]

tmp_all_param_ls <- initialize_galaina_param_and_input(config, fileConn = fileConn, Y = tmp_Y, Lambda = tmp_Lambda)

out_ls <- galaina_pipeline(tmp_all_param_ls)

# Compare
compare_v <- pcalg::compareGraphs(gl = out_ls[['out_graph']], gt = multi_out_ls[[i_run]][['original_graph']])
# compare_v['tpr']
# compare_v['fpr']
# compare_v['tdr']
cat(c(compare_v, "\n"))
shd_n <- pcalg::shd(out_ls[['out_graph']], multi_out_ls[[i_run]][['original_graph']])
cat(c(shd_n, "\n"))

# Warning messages:
# 1: In skeleton(suffStat, indepTest, alpha, labels = labels, method = skel.method,  :
# NAs introduced by coercion to integer range

    
# }
# if (sys.nframe() == 0) {
#   main()
# }
