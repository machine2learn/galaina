# Title     : Testing Pipeline for Galaina engine
# Objective : Generate data or load data from file; discover causal graphs; compare those graphs with ground truth graph
# Created by: fabio
# Created on: 02/03/2021

# rm(list = ls())

library(pcalg)
library(infotheo)
library(BDgraph)
library(sbgcop)

library(ConfigParser)

library(rlist)
# from demo_copulaFactorPC.R

source(file.path(getwd(), 'for_testing/artificial_generator_data_and_factor.R'))

source(file.path(getwd(), 'R_code/backend_galaina_final_part.R'))

# if (sys.nframe() == 0) {
main <- function(
  path_to_artificial_data = file.path(getwd(), 'artificial_data'),
  generator_parameters = NULL,
  overwrite_galaina_output = TRUE,
  measure_performance = TRUE,
  path_to_config = file.path(getwd(), 'user_data/test/AsiaSmall/config_5/config.ini'),
  tag_generator = NULL,
  tag_method = NULL
)
{
  # TODO - make a config.ini paying attention at:
  # input_path_directed_edges_blacklist

  cat(c("Initialize parameters", "\n"))
  fileConn <- NULL
  # generate_data_b <- is.null(path_to_artificial_data) | file_test( '-d', path_to_artificial_data)
  generate_data_b <- file_test( '-d', path_to_artificial_data)
  load_generated_data_b <- file_test( '-f', path_to_artificial_data)
  stopifnot(generate_data_b == !(load_generated_data_b))  # input must be path to file or folder
  # stopifnot(generate_data_b | load_generated_data_b)

  if (is.null(generator_parameters) & generate_data_b) {
    generator_parameters <- list(dataset_n = 3)
  }

  config <- ConfigParser$new(Sys.getenv(), optionxform = identity)
  config$read(filepath = path_to_config)

  if (generate_data_b) {
    cat(c(glue::glue("Generate {generator_parameters['dataset_n']} artificial datasets"), "\n"))
    multi_out_ls <- do.call('multi_artificial_generator_data_and_factor', generator_parameters)
    # multi_out_ls <- multi_artificial_generator_data_and_factor(dataset_n = dataset_n )
    name_artificial_data_output <- glue::glue("{multi_out_ls[[1]][['timestamp']]}.rdata")
    path_to_artificial_data_output <- file.path(path_to_artificial_data, name_artificial_data_output)
    print(path_to_artificial_data_output)
    list.save(x = multi_out_ls, file = path_to_artificial_data_output)
  } else {
    cat(c("Load artificial data", "\n"))
    path_to_artificial_data_output <- path_to_artificial_data
    multi_out_ls <- list.load(file = path_to_artificial_data_output)
  }

  # for (i_run in 1:length(multi_out_ls)) {
  #
  # }

  # if (gen)
  cat(c("Run Galaina", "\n"))
  load_copula_factor_functions(config)
  # Loop
  # i_run <- 1
  for (i_run in 1:length(multi_out_ls)) {
    tmp_Y <- multi_out_ls[[i_run]][['for_inference']][['Y']]
    tmp_Lambda <- multi_out_ls[[i_run]][['for_inference']][['Lambda']]

    tmp_all_param_ls <- initialize_galaina_param_and_input(config, fileConn = fileConn, Y = tmp_Y, Lambda = tmp_Lambda)
    if (overwrite_galaina_output){
      out_ls <- galaina_pipeline(tmp_all_param_ls)
    }

    # Compare
    compare_v <- pcalg::compareGraphs(gl = out_ls[['out_graph']], gt = multi_out_ls[[i_run]][['original_graph']])
    for (j_metric in names(compare_v)) {
      cat(c(glue::glue("{j_metric}: {compare_v[j_metric]}"), "\n"))
    }
    # compare_v['tpr']
    # compare_v['fpr']
    # compare_v['tdr']
    # cat(c(compare_v, "\n"))
    shd_n <- pcalg::shd(out_ls[['out_graph']], multi_out_ls[[i_run]][['original_graph']])
    # cat(c(shd_n, "\n"))

    cat(c(glue::glue("Structural Hamming Distance (SHD): {shd_n}"),  "\n"))
    # END - Loop

  }


  # Warning messages:
  # 1: In skeleton(suffStat, indepTest, alpha, labels = labels, method = skel.method,  :
  # NAs introduced by coercion to integer range

}

if (sys.nframe() == 0) {
  main()
}