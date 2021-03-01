#!/usr/bin/env Rscript

# This is the code called with the default config.ini

library(utils)
library(tools)  # to get file extension
library(dplyr)  # actually I load it only for filter
library(ggplot2)  # to save in SVG on R installed via homebrew  https://stackoverflow.com/questions/50767445/how-to-fix-failed-to-load-cairo-dll-in-r
library(ggplotify)  # to convert plot to grob objects
library(grid)
library(stringi)  # to check if strings are empty
library(sbgcop)
library(infotheo)
library(BDgraph)
library(pcalg)  # depends on bioconductor
library(bnlearn)
library(ConfigParser)  # GPL -
library(Rgraphviz)  # for using toDot
# library(configr)  # alternative to package above

# MAIN
if (sys.nframe() == 0) {
  messaging_b <- TRUE
  data_df <- NULL
  Lambda <- NULL
  # these are used for debugging (I don't need to pass Python stuff)
  #
  # args <- c("/Volumes/encrypted_data_sparse/data_for_demo/first_config.ini")
  # args <- c("/Volumes/encrypted_data_sparse/data_for_demo/first_unquoted_blacklist.ini")
  # args <- c("/Volumes/encrypted_data_sparse/data_for_demo/first_no_blacklist.ini")
  # After merging snap_other <- snap_teacher + snap_mother + snap_father
  # args <- c("/Users/fabio/projects/aggressotype/code/software/software_demo/config_estonian_restricted_demo.ini")
  # args <- c("/Users/fabio/projects/aggressotype/code/software/software_demo/config_estonian_restricted_demo_less_null.ini")
  #
  # args <- c("/Users/fabio/projects/aggressotype/code/software/software_demo/config_estonian_restricted_demo_relaxed_less_null.ini")

  if (messaging_b) {
    fileConn <- file("output.txt", 'a')  # TODO get it from config
    sink(file = fileConn, append = TRUE, type = 'message')  # future: unset type and use sink for managing output. Modify print_and_append_to_log because it could become redundant
  } else {
    fileConn <- NULL
  }

  # Read input config file and get all the arguments from that
  args <- commandArgs(trailingOnly = TRUE)  # taking config path argument from Python
  config <- ConfigParser$new(Sys.getenv(), optionxform = identity)
  config$read(filepath = args[1])
  # IDEA: make this piece of code (till loop) a function:
  # Input: (config_obj, messaging_b, data_df = NULL, Lambda = NULL)
  # Processing:
  # Output: pc_parameters_ls, infer_copula_param_ls, all_infer_param_ls, factor_names
  # BUT we need a special script for testing to specify input

  # args <- c("/Volumes/encrypted_data_sparse/demo_output/test_background_1000_bootstrap_3_merged_data.csv",
  #           "/Volumes/encrypted_data_sparse/demo_output/test_background_1000_bootstrap_3_merged_factor_loading_matrix.csv",
  #           "/Volumes/encrypted_data_sparse/demo_output/test_background_1000_bootstrap_1_causal_discovery.pdf",
  #           "1000", "500", "0",
  #           # "/Volumes/encrypted_data_sparse/demo_output/test_background_1000_bootstrap_3_correlation_matrix.csv",
  #           "/Volumes/encrypted_data_sparse/demo_output/test_background_1000_bootstrap_1_PcAlgo.rds",
  #           "/Volumes/encrypted_data_sparse/data_for_demo/background_knoledge.csv",
  #           "0"
  #           # "70"  # we need to store here the sample size
  #           )
  #           # "/Volumes/encrypted_data_sparse/data_for_demo/background_knoledge_exp.csv",  # edited edge
  column_separator_str <- config$get(option = "column_separator_str", section = "input_paths_and_related_parameters")  #',' <- "NULL"
  path_r_infer_copula_factor_script <- config$get(option = "path_r_infer_copula_factor_script", section = "r_front_end")  #',' <- "NULL"

  source(path_r_infer_copula_factor_script)
  # source('/Users/fabio/projects/aggressotype/code/software/software_demo/R_code/my_inferCopulaFactorModel.R')
  # source('my_inferCopulaFactorModel.R')

  #### 1. Load input ####
  # Load data_matrix, factor_loading_matrix
  # data_df <- read.csv(file = '/Volumes/encrypted_preprocessed_data/data_for_demo/intersection/data_df.csv', row.names = 1, na.strings=c("",".","NA"))
  # factor_loading_df <- read.csv(file = '/Volumes/encrypted_preprocessed_data/data_for_demo/intersection/factor_loading_df.csv', row.names = 1,  na.strings=c("",".","NA"))

  if (is.null(data_df)) {
    data_df <- read.csv(
      file = config$get(option = "output_path_merged_data", section = "output_paths"),
      sep = column_separator_str, row.names = 1, na.strings = c("", ".", "NA")
    )
  }

  if (is.null(Lambda)) {
  # if (is.null(factor_loading_df)) {
    factor_loading_df <- read.csv(
      # file = args[2],
      file = config$get(option = "output_path_merged_factor_model_loading", section = "output_paths"),
      sep = column_separator_str, row.names = 1, na.strings = c("", ".", "NA")
    )
    Lambda <- data.matrix(factor_loading_df)
    factor_names <- colnames(factor_loading_df)
  } else {
    factor_names <- 1:ncol(Lambda)
  }

  # load(file = '/Volumes/encrypted_preprocessed_data/demo_aggressotype_transfer.rda')
  # data_matrix <- as.numeric(data_matrix)
  # factor_loading_matrix <- as.numeric(factor_loading_matrix)

  # print(typeof(data_matrix))
  # print(typeof(factor_loading_matrix))

  gibbs_sampling_n <- as.integer(
    config$getfloat(option = "gibbs_sampling_n", section = "copula_factor_algorithm")
  )  # in INI e.g. 1000
  gibbs_burn_in_n <- as.integer(
    config$getfloat(option = "gibbs_burn_in_n", section = "copula_factor_algorithm")
  )  # in INI e.g. 500
  #
  gibbs_first_random_seed_n <- as.integer(
    config$getfloat(option = "gibbs_first_random_seed_n", section = "copula_factor_algorithm")
  )     # 365
  gibbs_random_seed_update_parameter_n <- as.integer(
    config$getfloat(option = "gibbs_random_seed_update_parameter_n", section = "copula_factor_algorithm")
  )

  original_bootstrap_n <- as.integer(
    config$getfloat(option = "bootstrap_n", section = "edge_weight_algorithm")
    # args[6]
  )  # in INI e.g. 500
  causal_discovery_algorithm_run_n <- max(original_bootstrap_n, 1)  # if original_bootstrap_n ==0, we know that no bootstrap must be performed

  # perform_bootstrap_b <- (original_bootstrap_n >= 1)
  no_perform_bootstrap_b <- (original_bootstrap_n == 0)
  perform_bootstrap_b <- !no_perform_bootstrap_b
  # if (perform_bootstrap_b) {
  bootstrap_first_random_seed_n <- as.integer(
    config$getfloat(option = "bootstrap_first_random_seed_n", section = "edge_weight_algorithm")
    # args[6]
  ) #100
  bootstrap_random_seed_update_parameter_n <- as.integer(
    config$getfloat(option = "bootstrap_random_seed_update_parameter_n", section = "edge_weight_algorithm")
    # args[6]
  )
  bootstrap_random_seed_n <- bootstrap_first_random_seed_n
  # }
  input_path_directed_edges_blacklist <- config$get(option = "input_path_directed_edges_blacklist", fallback = "", section = "input_paths_and_related_parameters")   # args[8]

  print_and_append_to_log(c("Gibbs total samples:", gibbs_sampling_n, "\n"), fileConn)
  print_and_append_to_log(c("Gibbs burn-in samples:", gibbs_burn_in_n, "\n"), fileConn)

  # TODO future: use whitelist and blacklist
  # if (file.exists(input_path_directed_edges_blacklist)) {
  if (!stri_isempty(input_path_directed_edges_blacklist)) {
    dir_arc_blacklist_mat <- as.matrix(read.csv(file = input_path_directed_edges_blacklist, sep = column_separator_str, row.names = NULL))
    # input_bgk_whitelist_from_to_mat <- as.matrix(read.csv(file = args[9], sep = column_separator_str, row.names = NULL))
  }
  # if (file.exists(args[10])) {
  #   input_bgk_blacklist_from_to_mat <- as.matrix(read.csv(file = args[10], sep = column_separator_str, row.names = NULL))
  # }
  # input_bgk_from_to_mat <- cbind(y, x)  # better use form above because this is to from and now we useded the ones with rev columns at the same time
  ##

  # causal_discovery_observation_n <- as.integer(config$getfloat(option = "causal_discovery_observation_n", section = "pc_algorithm"))

  # if (original_bootstrap_n >= 1) {
  #   bootstrap_random_seed_n <- bootstrap_first_random_seed_n
  # }
  #### 2. Perform bootstrap, inference, causal discovery ####
  ## some checks
  stopifnot(gibbs_burn_in_n < gibbs_sampling_n)  # burn-in must be lower
  stopifnot(all(sapply(data_df, is.numeric)))  # check all columns in data df are numeric
  # stopifnot(all(sapply(factor_loading_df, is.numeric)))  # check all columns in factor df are numeric
  stopifnot(all(sapply(Lambda, is.numeric)))  # TODO see if it works - check all columns in factor df are numeric

  run2pcalgo_ls <- vector(mode = "list", length = causal_discovery_algorithm_run_n)
  run2pcalgo_bad_ls <- run2pcalgo_ls  # Stores the ones we removed
  run2suffstat_ls <- run2pcalgo_ls

  # bootstrapped_graphnel_list <- vector(mode = "list", length = causal_discovery_algorithm_run_n)

  # TODO if need this parameter in the future, fixedGaps and fixedEdges will be logical matrices
  pc_fixedGaps <- config$get(option = "fixedGaps", section = "pc_algorithm")
  if (identical(pc_fixedGaps, "NULL")) {
    pc_fixedGaps <- NULL
  }
  pc_fixedEdges <- config$get(option = "fixedEdges", section = "pc_algorithm")
  if (identical(pc_fixedEdges, "NULL")) {
    pc_fixedEdges <- NULL
  }

  # cat(config$data$pc_algorithm$indepTest, "\n")
  indep_test_parameter_str <- config$get(option = "indepTest", section = "pc_algorithm")
  output_path_suffstat <- config$get(option = "output_path_suffstat", fallback = "", section = "output_paths")

  output_path_pc_algo_obj <- config$get(option = "output_path_pc_algo_obj", fallback = "", section = "output_paths")

  pc_parameters_ls <- list(
    suffStat = list(),
    indepTest = get(indep_test_parameter_str),  # now it pastes here the code of the function
    # indep_test_parameter_str, #<- indep_test_parameter_str[[1]], # substitute(indep_test_parameter_str),
    # indepTest = config$get(option = "indepTest", section = "pc_algorithm"),
    alpha = config$getfloat(option = "alpha", section = "pc_algorithm"),
    labels = factor_names,
    numCores = as.integer(config$getfloat(option = "numCores", section = "pc_algorithm")),
    verbose = config$getboolean(option = "verbose", section = "pc_algorithm"),
    fixedGaps = pc_fixedGaps,
    fixedEdges = pc_fixedEdges,
    NAdelete = config$getboolean(option = "NAdelete", section = "pc_algorithm"),
    m.max = config$get(option = "m.max", section = "pc_algorithm"),
    u2pd = config$get(option = "u2pd", section = "pc_algorithm"),
    # u2pd = "relaxed",  # impled by solve.confl = TRUE
    skel.method = config$get(option = "skel.method", section = "pc_algorithm"),
    conservative = config$getboolean(option = "conservative", section = "pc_algorithm"),
    # conservative = TRUE,
    # Choose either conservative PC or majority rule PC!
    maj.rule = config$getboolean(option = "maj.rule", section = "pc_algorithm"),
    # maj.rule = TRUE,
    solve.confl = config$getboolean(option = "solve.confl", section = "pc_algorithm")
    # solve.confl = TRUE
  )

  # TODO maybe in the future we require to have at least causal_discovery_algorithm_run_n algorithm, so we replace this for with a while (iRun_n <=  causal_discovery_algorithm_run_n)
  # This loop outputs: run2suffstat_ls, run2pcalgo_ls, run2pcalgo_bad_ls, bootstrap_random_seed_n
  # We could turn the loop core into a function taking as input
  # run2suffstat_ls, run2pcalgo_ls, run2pcalgo_bad_ls, gibbs_sampling_n, factor_loading_df
  # AND causal_discovery_algorithm_run_n
  # AND
  # no_perform_bootstrap_b, fileConn, data_df, bootstrap_random_seed_n, bootstrap_random_seed_update_parameter_n
  # AND
  # gibbs_sampling_n, odens_n, first_random_seed, random_seed_update_parameter,
  # AND other stuff adn outputing I think those 3 and
  # original_bootstrap_n, fileConn, bootstrap_random_seed_n

  odens_n <- 1  # save only every odens_n
  #

  infer_copula_param_ls <- list(
    Y = list(NULL),
    Lambda = Lambda,  # data.matrix(factor_loading_df),
    nsamp = gibbs_sampling_n,
    odens = odens_n,  # Store each Gibbs sampling output
    first_random_seed = gibbs_first_random_seed_n,
    random_seed_update_parameter = gibbs_random_seed_update_parameter_n,
    output_intermediate = output_path_suffstat,
    run_id = list(NULL),  # save intermediate Sigma.psamp for each bootstrap/run
    fileConn = list(NULL)
  )
  if (messaging_b) {
    infer_copula_param_ls[['fileConn']] <- fileConn
  }

  # factor_n <- length(factor_names)
  # throw_away_odens_n <- floor(gibbs_burn_in_n / odens_n)

  all_infer_param_ls <- list(
    data_df = data_df,
    # factor_n = length(factor_names),
    factor_n = length(pc_parameters_ls[['labels']]),
    throw_away_odens_n = floor(gibbs_burn_in_n / odens_n),
    infer_copula_param_ls = infer_copula_param_ls,
    pc_parameters_ls = pc_parameters_ls,
    # output_path_suffstat = output_path_suffstat,
    output_path_pc_algo_obj = output_path_pc_algo_obj,
    run2pcalgo_ls = run2pcalgo_ls,
    run2pcalgo_bad_ls = run2pcalgo_bad_ls,
    run2suffstat_ls = run2suffstat_ls,
    causal_discovery_observation_n = as.integer(
      config$getfloat(option = "causal_discovery_observation_n", section = "pc_algorithm")
    ),
    # causal_discovery_observation_n,
    iRun_n = list(NULL),
    causal_discovery_algorithm_run_n = causal_discovery_algorithm_run_n,
    perform_bootstrap_b = perform_bootstrap_b,
    bootstrap_random_seed_n = bootstrap_random_seed_n,
    bootstrap_random_seed_update_parameter_n = bootstrap_random_seed_update_parameter_n,
    fileConn = infer_copula_param_ls[['fileConn']]
  )
  # path_to_bn_strength_obj <- config$get(
  #   option = "output_path_bn_strength_obj", fallback = "", section = "output_paths"
  # )
  # path_to_avg_bn_obj <- config$get(
  #   option = "output_path_avg_bn_obj", fallback = "", section = "output_paths"
  # )
  # main_plot_title_str <- config$get(
  #   option = "core_plot_title_str", fallback = "", section = "plot_and_display"
  # )
  # tmp_fig <- config$get(option = "output_path_fig", section = "output_paths")

  post_bootstrap_processing_param_ls <- list(
    # path_to_bn_strength_obj = path_to_bn_strength_obj,
    # path_to_avg_bn_obj = path_to_avg_bn_obj,
    # main_plot_title_str = main_plot_title_str,
    # tmp_fig = tmp_fig,
    path_to_bn_strength_obj = config$get(
      option = "output_path_bn_strength_obj", fallback = "", section = "output_paths"
    ),
    path_to_avg_bn_obj = config$get(
      option = "output_path_avg_bn_obj", fallback = "", section = "output_paths"
    ),
    main_plot_title_str = config$get(
      option = "core_plot_title_str", fallback = "", section = "plot_and_display"
    ),
    tmp_fig = config$get(option = "output_path_fig", section = "output_paths"),
    factor_names = pc_parameters_ls[['labels']],  # factor_names  # also in pc_param_ls
    perform_bootstrap_b = all_infer_param_ls[['perform_bootstrap_b']]
  )
  # return(pc_parameters_ls, infer_copula_param_ls, all_infer_param_ls, post_bootstrap_processing_param_ls, output_path_suffstat)
  #
  # )
  #
  tmp_all_infer_param_ls <- all_infer_param_ls
  for (iRun_n in 1:causal_discovery_algorithm_run_n) {
    # tmp_out_ls <- infer_covariance_and_graph(
    #   data_df, factor_n, throw_away_odens_n,
    #   infer_copula_param_ls,
    #   pc_parameters_ls,
    #   output_path_pc_algo_obj,
    #   run2pcalgo_ls,
    #   run2pcalgo_bad_ls,
    #   run2suffstat_ls,
    #   causal_discovery_observation_n,
    #   iRun_n,
    #   causal_discovery_algorithm_run_n,
    #   perform_bootstrap_b,
    #   bootstrap_random_seed_n,
    #   bootstrap_random_seed_update_parameter_n,
    #   fileConn
    # )
    # bootstrap_random_seed_n <- tmp_out_ls[['bootstrap_random_seed_n']]
    # run2pcalgo_ls <- tmp_out_ls[['run2pcalgo_ls']]
    # run2pcalgo_bad_ls <- tmp_out_ls[['run2pcalgo_bad_ls']]
    # run2suffstat_ls <- tmp_out_ls[['run2suffstat_ls']]
    #
    tmp_all_infer_param_ls[['iRun_n']] <- iRun_n
    tmp_out_ls <- do.call("infer_covariance_and_graph", tmp_all_infer_param_ls)
    tmp_all_infer_param_ls[['bootstrap_random_seed_n']] <- tmp_out_ls[['bootstrap_random_seed_n']]
    tmp_all_infer_param_ls[['run2pcalgo_ls']] <- tmp_out_ls[['run2pcalgo_ls']]
    tmp_all_infer_param_ls[['run2pcalgo_bad_ls']] <- tmp_out_ls[['run2pcalgo_bad_ls']]
    tmp_all_infer_param_ls[['run2suffstat_ls']] <- tmp_out_ls[['run2pcalgo_ls']]
  }
  run2pcalgo_ls <- tmp_out_ls[['run2pcalgo_ls']]
  run2pcalgo_bad_ls <- tmp_out_ls[['run2pcalgo_bad_ls']]
  run2suffstat_ls <- tmp_out_ls[['run2suffstat_ls']]

  output_path_pc_algo_obj <- all_infer_param_ls[['output_path_pc_algo_obj']]
  # output_path_suffstat <- all_infer_param_ls[['output_path_suffstat']]

  if (!stri_isempty(output_path_pc_algo_obj)) {
    saveRDS(run2pcalgo_ls, file = output_path_pc_algo_obj)
  }
  if (!stri_isempty(output_path_suffstat)) {
    saveRDS(run2suffstat_ls, file = output_path_suffstat)
  }
  # return(list(
  #   run2pcalgo_ls = run2pcalgo_ls,
  #   run2pcalgo_bad_ls = run2pcalgo_bad_ls,
  #   run2suffstat_ls = run2suffstat_ls
  # ))
  cat('\n')

  # Save stuff generated in the loop above
  out_graph <- merge_discovered_graphs(post_bootstrap_processing_param_ls, run2pcalgo_ls)

  dev.off()
  print_and_append_to_log(c("Completed", "\n"), fileConn)
  if (messaging_b) {
    sink(type = "message")  # reset message and output sink
    close(fileConn)
  }
}