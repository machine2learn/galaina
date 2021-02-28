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
library(ConfigParser)  # GPL-
library(Rgraphviz)  # for using toDot
# library(configr)  # alternative to package above

# MAIN
if (sys.nframe() == 0) {
  messaging_b <- TRUE
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
  }

  # Read input config file and get all the arguments from that
  args <- commandArgs(trailingOnly = TRUE)  # taking config path argument from Python
  config <- ConfigParser$new(Sys.getenv(), optionxform = identity)
  config$read(filepath = args[1])
  # IDEA: make this piece of code (till loop) a function:
  # Input: (config file, messaging_b, optional data_df, optional factor_loading_df)
  # Processing: initialize stuff
  # Output: pc_parameters_ls, infer_copula_param_ls, bootstrap_param_ls
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

  data_df <- read.csv(
    file = config$get(option = "output_path_merged_data", section = "output_paths"),
    sep = column_separator_str, row.names = 1, na.strings = c("", ".", "NA")
  )

  factor_loading_df <- read.csv(
    # file = args[2],
    file = config$get(option = "output_path_merged_factor_model_loading", section = "output_paths"),
    sep = column_separator_str, row.names = 1, na.strings = c("", ".", "NA")
  )

  factor_names <- colnames(factor_loading_df)
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
  causal_discovery_algorithm_run_n <- max(original_bootstrap_n, 1)  # if original_bootstrap_n ==0, we know that no boostrap must be performed

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

  causal_discovery_observation_n <- as.integer(config$getfloat(option = "causal_discovery_observation_n", section = "pc_algorithm"))

  # if (original_bootstrap_n >= 1) {
  #   bootstrap_random_seed_n <- bootstrap_first_random_seed_n
  # }
  #### 2. Perform bootstrap, inference, causal discovery ####
  ## some checks
  stopifnot(gibbs_burn_in_n < gibbs_sampling_n)  # burn-in must be lower
  stopifnot(all(sapply(data_df, is.numeric)))  # check all columns in data df are numeric
  stopifnot(all(sapply(factor_loading_df, is.numeric)))  # check all columns in factor df are numeric

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
    Lambda = data.matrix(factor_loading_df),
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

  factor_n <- length(factor_names)
  throw_away_odens_n <- floor(gibbs_burn_in_n / odens_n)

  # bootstrap_iterarion_param_ls <- list(
  #
  # )
  for (iRun_n in 1:causal_discovery_algorithm_run_n) {
    tmp_out_ls <- infer_covariance_and_graph(
      data_df, factor_n, throw_away_odens_n,
      infer_copula_param_ls,
      pc_parameters_ls,
      output_path_pc_algo_obj,
      run2pcalgo_ls,
      run2pcalgo_bad_ls,
      run2suffstat_ls,
      causal_discovery_observation_n,
      iRun_n,
      causal_discovery_algorithm_run_n,
      perform_bootstrap_b,
      bootstrap_random_seed_n,
      bootstrap_random_seed_update_parameter_n,
      fileConn
    )
    bootstrap_random_seed_n <- tmp_out_ls[['bootstrap_random_seed_n']]
    run2pcalgo_ls <- tmp_out_ls[['run2pcalgo_ls']]
    run2pcalgo_bad_ls <- tmp_out_ls[['run2pcalgo_bad_ls']]
    run2suffstat_ls <- tmp_out_ls[['run2suffstat_ls']]
  }
  run2pcalgo_ls <- tmp_out_ls[['run2pcalgo_ls']]
  run2pcalgo_bad_ls <- tmp_out_ls[['run2pcalgo_bad_ls']]
  run2suffstat_ls <- tmp_out_ls[['run2suffstat_ls']]
  #
  #   ## 2.1 Bootstrap
  #   # Input:
  #   # - iRun_n, causal_discovery_algorithm_run_n, no_perform_bootstrap_b, bootstrap_random_seed_n, bootstrap_random_seed_update_parameter_n,
  #   # - data_df,
  #   # - infer_copula_param_ls, pc_parameters_ls, identity
  #   # - throw_away_odens_n, factor_n, causal_discovery_observation_n,
  #   # - output_path_pc_algo_obj
  #   # - (fileConn),
  #   # Output:
  #   # - bootstrap_random_seed_n
  #   # return(list(run2suffstat_ls = run2suffstat_ls, run2pcalgo_ls = run2pcalgo_ls, run2pcalgo_bad_ls = run2pcalgo_bad_ls))
  #
  #   if (no_perform_bootstrap_b) {
  #     # cat("No bootstrap performed \n")
  #     print_and_append_to_log(c("No bootstrap performed \n"), fileConn)
  #     tmp_run_data_df <- data_df
  #   } else {
  #     set.seed(bootstrap_random_seed_n) # for testing set seed to
  #     print_and_append_to_log(c("Current bootstrap sample:", iRun_n, "of", causal_discovery_algorithm_run_n, "\n"), fileConn)
  #     print_and_append_to_log(c("Random seed bootstrap sample:", bootstrap_random_seed_n, "\n"), fileConn)
  #     iteration_rows_id <- sample(nrow(data_df), replace = TRUE)
  #     bootstrap_random_seed_n <- update_random_seed(bootstrap_random_seed_n, bootstrap_random_seed_update_parameter_n)  # FG update random seed
  #     tmp_run_data_df <- data_df[iteration_rows_id, ]
  #     # write.csv(
  #     #   tmp_run_data_df,
  #     #   file = '/Users/fabio/rpq_230.csv',
  #     #   row.names=FALSE,
  #     #   na=""
  #     # )
  #   }
  #
  #   ## 2.2 Perform inference - first seed of gibbs sampling is always the same, but data changes
  #   print_and_append_to_log("Perform inference\n", fileConn)
  #   tmp_infer_copula_param_ls <- infer_copula_param_ls
  #   tmp_infer_copula_param_ls[['Y']] <- data.matrix(tmp_run_data_df)
  #   tmp_infer_copula_param_ls[['run_id']] <- iRun_n
  #   if (no_perform_bootstrap_b) {  # Used for testing, but in theory the repeat would be enough
  #     tmp_cop_fac_obj <- do.call("my_inferCopulaFactorModel", tmp_infer_copula_param_ls)
  #   } else {
  #     repeat {
  #       tmp_cop_fac_obj <- tryCatch(
  #         tmp_cop_fac_obj <- do.call("my_inferCopulaFactorModel", tmp_infer_copula_param_ls),
  #         error = identity
  #       )
  #       if (!is(tmp_cop_fac_obj, "error") || no_perform_bootstrap_b) {
  #         break
  #       }
  #       if (is(tmp_cop_fac_obj, "error")) {
  #         print_and_append_to_log(c("Error:", tmp_cop_fac_obj[[1]], "\n"), fileConn)
  #         traceback()
  #         #print_and_append_to_log(c("Error", "\n"), fileConn)
  #         #for (j_el in tmp_cop_fac_obj) {
  #         #  if (is.list(j_el)) {
  #         #    print(j_el)
  #         #    #out_txt <- paste0(unlist(j_el), collapse = "\n")
  #         #    #for (k_el in j_el) {
  #         #    print_and_append_to_log(c(out_txt, "\n"), fileConn)
  #         #    #}
  #         #  } else {
  #         #    print_and_append_to_log(c(j_el, "\n"), fileConn)
  #         #  }
  #         #}
  #       }
  #       set.seed(bootstrap_random_seed_n) # for testing set seed to
  #       print_and_append_to_log(c("Repeat current bootstrap sample:", iRun_n, "of", causal_discovery_algorithm_run_n, "\n"), fileConn)
  #       print_and_append_to_log(c("Random seed bootstrap sample:", bootstrap_random_seed_n, "\n"), fileConn)
  #       iteration_rows_id <- sample(nrow(data_df), replace = TRUE)
  #       bootstrap_random_seed_n <- update_random_seed(bootstrap_random_seed_n, bootstrap_random_seed_update_parameter_n)  # FG update random seed
  #       tmp_run_data_df <- data_df[iteration_rows_id, ]
  #     }
  #   }
  #
  #
  #   # cat(dim(tmp_cop_fac_obj$Sigma.psamp), "\n")
  #   # Extract samples of the correlation matrix over latent variables, ignoring the first samples (burn-in)
  #   # C_samples <- tmp_cop_fac_obj$Sigma.psamp[1:factor_n, 1:factor_n, (gibbs_burn_in_n + 1) : gibbs_sampling_n]
  #   if (throw_away_odens_n < 1) {
  #     C_samples <- tmp_cop_fac_obj$Sigma.psamp[1:factor_n, 1:factor_n, ]
  #   } else {
  #     C_samples <- tmp_cop_fac_obj$Sigma.psamp[1:factor_n, 1:factor_n, -(1:throw_away_odens_n)]
  #   }
  #   # the posterior mean
  #   C <- apply(C_samples, c(1, 2), mean)
  #
  #   # Allows user to input a fix value to causal_discovery_observation_n
  #   tmp_causal_discovery_observation_n <- causal_discovery_observation_n
  #   if (causal_discovery_observation_n == 0) {
  #     C_sd <- apply(C_samples, c(1, 2), sd)
  #     # effective sample size
  #     C_ess <- ((1 - C^2) / C_sd)^2
  #     # effective_observation_n <- mean(C_ess[upper.tri(C_ess)])
  #     tmp_causal_discovery_observation_n <- mean(C_ess[upper.tri(C_ess)])
  #   }
  #
  #   ## export to CSV
  #   # write.csv(C, file = args[7])
  #
  #   # set.seed(random_seed_n)
  #   #
  #   #### 2.3 Causal discovery ####
  #   # MEMO we could directly store the bnlearn version
  #   ## call the order independent version of the standard PC algorithm
  #
  #   tmp_suffstat_ls <- list(C = C, n = tmp_causal_discovery_observation_n)
  #   run2suffstat_ls[[iRun_n]] <- tmp_suffstat_ls
  #   # For parallelization, we could just make a function returning tmp_suffstat_ls and store it into run2suffstat_ls.
  #   # After that, in another loop, we can run the causal discovery algorithm and store
  #   # - run2pcalgo_ls
  #   # - run2pcalgo_bad_ls
  #
  #   tmp_pc_parameters_ls <- pc_parameters_ls
  #   tmp_pc_parameters_ls[['suffStat']] <- tmp_suffstat_ls
  #   # stopifnot(!identical(old_pc_parameters_ls$suffStat, pc_parameters_ls$suffStat))
  #   tmp_graph_cfpc <- do.call("pc", tmp_pc_parameters_ls)
  #
  #   # WARNING we might have to turn off the edge assignment part
  #
  #   # DEBUG begin
  #   # Restrict to arc that can be enforced on this graph: they must correspond to egdes that are undirected
  #   # Check if edge is there
  #   tmp_blacklist_arcs_absent_b <- TRUE
  #   if (file.exists(input_path_directed_edges_blacklist)) {
  #     tmp_amat_cpdag <- as(tmp_graph_cfpc, "amat")
  #     # DAG/CPDAG (format "cpdag"). Directed egde {from} --> {to} is present <=> amat[{from}, {to}] == 0 AND amat[{to}, {from}] == 1
  #     tmp_cpdag_blacklist_arcs_present_bv <- (tmp_amat_cpdag[dir_arc_blacklist_mat] == 0) & (tmp_amat_cpdag[dir_arc_blacklist_mat[, rev(colnames(dir_arc_blacklist_mat))]] == 1)
  #     # MAG/DAG (format "pag"). Directed egde {from} --> {to} is present <=> amat[{from}, {to}] == 2 AND amat[{to}, {from}] == 3
  #     tmp_pag_blacklist_arcs_present_bv <- (tmp_amat_cpdag[dir_arc_blacklist_mat] == 2) & (tmp_amat_cpdag[dir_arc_blacklist_mat[, rev(colnames(dir_arc_blacklist_mat))]] == 3)
  #     tmp_blacklist_arcs_present_bv <- (tmp_cpdag_blacklist_arcs_present_bv | tmp_pag_blacklist_arcs_present_bv)
  #     # tmp_blacklist_arcs_present_bv <- (tmp_amat_cpdag[dir_arc_blacklist_mat] == 1) & (tmp_amat_cpdag[dir_arc_blacklist_mat[, rev(colnames(dir_arc_blacklist_mat))]] == 0)
  #     tmp_blacklist_arcs_absent_b <- !any(tmp_blacklist_arcs_present_bv)
  #   }
  #   # if (file.exists(args[9])) {
  #   #   tmp_arc_enforceable_bv <- as.logical(tmp_amat_cpdag[input_bgk_from_to_mat]) & as.logical(tmp_amat_cpdag[input_bgk_from_to_mat[, rev(colnames(input_bgk_from_to_mat))]])
  #   #   # tmp_arc_enforceable_bv <- as.logical(tmp_amat_cpdag[input_bgk_from_to_mat])
  #   #   # cat(tmp_arc_enforceable_bv, "\n")
  #   #   # If some enforceable arc are still there, enforce them
  #   #   if (any(tmp_arc_enforceable_bv)) {
  #   #     # Maybe check x and y are not empty
  #   #     with_bgk_amat_cpdag <- addBgKnowledge(
  #   #       gInput = tmp_amat_cpdag,
  #   #       # x = x[tmp_arc_enforceable_bv],
  #   #       # y = y[tmp_arc_enforceable_bv],
  #   #       x = input_bgk_from_to_mat[tmp_arc_enforceable_bv, "from"],
  #   #       y = input_bgk_from_to_mat[tmp_arc_enforceable_bv, "to"],
  #   #       verbose = TRUE,
  #   #       checkInput = FALSE  # because it was giving error
  #   #     )
  #   #     # cat(with_bgk_amat_cpdag, "\n")
  #   #     # cat(typeof(with_bgk_amat_cpdag), "\n")
  #   #   }
  #   # }
  #   # bootstrapped_graphnel_list[[iRun_n]] <- as(t(as(with_bgk_amat_cpdag, "matrix")), "graphNEL")  # strange we need to specify matrix
  #
  #   # In theory we would need to convert to graphNEL only if we add background knowledge, because as.bn would work with both
  #   # bootstrapped_graphnel_list[[iRun_n]] <- as(t(with_bgk_amat_cpdag), "graphNEL")
  #   if (tmp_blacklist_arcs_absent_b) {
  #     run2pcalgo_ls[[iRun_n]] <- tmp_graph_cfpc
  #   } else {
  #     run2pcalgo_bad_ls[[iRun_n]] <- tmp_graph_cfpc
  #     print_and_append_to_log(
  #       c("Current graph removed because it contained", sum(tmp_blacklist_arcs_present_bv), "blacklisted directed edges", "\n"),
  #       fileConn
  #     )
  #     stopifnot(sum(tmp_blacklist_arcs_present_bv) > 0)
  #     print_and_append_to_log(c(dir_arc_blacklist_mat[tmp_blacklist_arcs_present_bv]), fileConn)
  #   }
  #   # saveRDS(graph_cfpc, file = args[8])
  #   cat('\n\n')
  #   # Save temporary output
  #   if (!stri_isempty(output_path_pc_algo_obj)) {
  #     tmp_out_path <- paste0(
  #       file_path_sans_ext(output_path_pc_algo_obj), '_', iRun_n, '.', file_ext(output_path_pc_algo_obj)
  #     )
  #     saveRDS(run2pcalgo_ls, file = tmp_out_path)
  #   }
  # }
  # return(list(run2suffstat_ls = run2suffstat_ls, run2pcalgo_ls = run2pcalgo_ls, run2pcalgo_bad_ls = run2pcalgo_bad_ls))

  cat('\n')

  # Save stuff generated in the loop above
  # maybe TODO adjust it when no bootstap is done
  output_path_suffstat <- config$get(option = "output_path_suffstat", fallback = "", section = "output_paths")
  if (!stri_isempty(output_path_suffstat)) {
    saveRDS(run2suffstat_ls, file = output_path_suffstat)
  }
  output_path_pc_algo_obj <- config$get(option = "output_path_pc_algo_obj", fallback = "", section = "output_paths")
  if (!stri_isempty(output_path_pc_algo_obj)) {
    saveRDS(run2pcalgo_ls, file = output_path_pc_algo_obj)
  }

  if (no_perform_bootstrap_b) {
    bn_obj <- as.bn(run2pcalgo_ls[[1]], check.cycles = FALSE)  # We are ok with cycles
  } else {
    removed_graph_due_to_blacklist_bv <- sapply(X = run2pcalgo_ls, FUN = is.null)
    # cat("Graphs removed due to incompatibility with blacklist", sum(removed_graph_due_to_blacklist_bv), "\n")
    print_and_append_to_log(c("Graphs removed due to incompatibility with blacklist:", sum(removed_graph_due_to_blacklist_bv), "\n"), fileConn)
    bootstrapped_bnlearn_list <- lapply(
      X = run2pcalgo_ls[!removed_graph_due_to_blacklist_bv],  # Ignore NULL elements of the list
      # X = run2pcalgo_ls[!sapply(run2pcalgo_ls, is.null)],  # Ignore NULL elements of the list
      # X = run2pcalgo_ls,
      FUN = as.bn,
      check.cycles = FALSE   # We are ok with cycles
    )
    # saveRDS(bootstrapped_graphnel_list, file = args[8])
    # bootstrapped_bnlearn_list <- lapply(X = bootstrapped_graphnel_list, FUN = as.bn)
    bn_strength_obj <- custom.strength(
      networks = bootstrapped_bnlearn_list,
      nodes = factor_names,
      cpdag = FALSE  # If TRUE the (PDAG of) the equivalence class is used instead of the network structure itself. It should make it easier to identify score-equivalent arcs.
    )
    avg_bn_obj <- averaged.network(strength = bn_strength_obj, nodes = factor_names)  # it automatically gets the thresold from bn_strength_obj

    path_to_bn_strength_obj <- config$get(
      option = "output_path_bn_strength_obj", fallback = "", section = "output_paths"
    )
    if (!stri_isempty(path_to_bn_strength_obj)) {
      saveRDS(bn_strength_obj, file = path_to_bn_strength_obj)
      # Save to CSV files (thresholded and not thresholded)
      bn_strength_obj_df <- readRDS(path_to_bn_strength_obj)
      write.csv(
        bn_strength_obj_df,
        file = sub('.rds', '.csv', path_to_bn_strength_obj),
        row.names = FALSE,
        na = ""
      )
      write.csv(
        bn_strength_obj_df[bn_strength_obj_df$strength >= attr(bn_strength_obj_df, "threshold"),],
        file = sub('.rds', '_thresholded.csv', path_to_bn_strength_obj),
        row.names = FALSE,
        na = ""
      )
    }

    path_to_avg_bn_obj <- config$get(
      option = "output_path_avg_bn_obj", fallback = "", section = "output_paths"
    )
    if (!stri_isempty(path_to_avg_bn_obj)) {
      saveRDS(avg_bn_obj, file = path_to_avg_bn_obj)
    }
  }

  main_plot_title_str <- config$get(option = "core_plot_title_str", fallback = "", section = "plot_and_display")
  # show_graph_before_bgk_application_b <- config$getboolean(option = "show_graph_before_bgk_application_b", section = "plot_and_display")  # as.logical("False")  # TODO put in INI
  # show_graph_after_bgk_application_b <- config$getboolean(option = "show_graph_after_bgk_application_b", section = "plot_and_display")  # as.logical("True")  # TODO put in INI
  # show_graph_with_dashed_removed_b <- config$getboolean(option = "show_graph_with_dashed_removed_b", section = "plot_and_display") # as.logical("False")  # TODO put in INI

  # Check if there are some blacklist directed edges to be removed
  # if (original_bootstrap_n >= 1){
  #   dir_arcs_mat <- directed.arcs(avg_bn_obj)
  # } else {
  #   dir_arcs_mat <- directed.arcs(bn_obj)
  # }
  # dir_arcs_to_be_removed_n <- 0
  # if ((nrow(dir_arcs_mat) > 0) & (file.exists(input_path_directed_edges_blacklist)) & (nrow(dir_arc_blacklist_mat) > 0)) {  # if at least one od those mat is empty, there are no arcs to be removed
  #   dir_arcs_and_blacklist_mat <- rbind(unique(dir_arcs_mat), unique(dir_arc_blacklist_mat))
  #   dir_arcs_in_blacklist_mat <- dir_arcs_and_blacklist_mat[duplicated(dir_arcs_and_blacklist_mat), , drop = FALSE]
  #   dir_arcs_to_be_removed_n <- nrow(dir_arcs_in_blacklist_mat)
  # }
  #
  # make_after_removal_plot_b <- (show_graph_after_bgk_application_b & (dir_arcs_to_be_removed_n > 0))
  # if (!(make_after_removal_plot_b | show_graph_before_bgk_application_b)) {
  #   show_graph_before_bgk_application_b <- TRUE  # if there is no second plot to make, force to make the first one
  # }

  graphics.off()

  tmp_fig <- config$get(option = "output_path_fig", section = "output_paths")
  # pdf(
  #   config$get(option = "output_path_fig_pdf", section = "output_paths"),
  #   width=6, height=3, onefile = FALSE
  # )
  do.call(
    file_ext(tmp_fig),
    list(tmp_fig, width = 6, height = 3, onefile = FALSE)
  )

  # Used when we were enforcing balacklist at the end
  # par(mfrow=c(show_graph_before_bgk_application_b + make_after_removal_plot_b, 1))
  # par(mar=c(1,1,1,1))
  par(mfrow = c(1, 1))

  first_title <- main_plot_title_str  # 'Causal graph with all edges'
  # Used when we were enforcing balacklist at the end
  # if (make_after_removal_plot_b) {
  #   first_title <- paste(main_plot_title_str, "- including blacklist")  #'Before removing blacklist edges'
  # }
  # if (!(make_after_removal_plot_b) & (dir_arcs_to_be_removed_n == 0)) {
  #   first_title <- paste(main_plot_title_str, "- no blacklist effect")  # No edges were removed
  # }
  # second_title <- main_plot_title_str
  # if (make_after_removal_plot_b & show_graph_with_dashed_removed_b) {
  #   second_title <- paste(main_plot_title_str, "- blacklist dashed")  #'Before removing blacklist edges'
  # } else if (make_after_removal_plot_b & !show_graph_with_dashed_removed_b) {
  #   second_title <- paste(main_plot_title_str, "- excluding blacklist")  #'Before removing blacklist edges'
  # }

  # cat(first_title, "\n")
  # if (show_graph_before_bgk_application_b) {
  if (no_perform_bootstrap_b) {
    # TODO - output to var, convert and save it
    plot(run2pcalgo_ls[[1]], main = first_title)  # if
  } else {
    # TODO - output to var, convert to bn <- strength (if needed), convert to DOT format or other and save to path specified in INI
    graph_nel <- strength.plot(x = avg_bn_obj, strength = bn_strength_obj, main = first_title)
        # threshold = attributes(bn_strength_obj)$threshold, # No need to specify the thershold
    # fix edge
    # Aug 2019 -
    # strength.plot() now saves arc strengths as weights in the graphNEL object it returns.
    tmp_key_v <- names(graph_nel@renderInfo@edges[["lwd"]])
    lwd_key2data_key <- as.list(setNames(object = gsub("[~]", "|", tmp_key_v), nm = tmp_key_v))
    for (iName in names(lwd_key2data_key)) {
      # graph_nel@edgeData@data[[lwd_key2data_key[[iName]]]][["penwidth"]] <- graph_nel@edgeData@data[iName]]
      graph_nel@edgeData@data[[lwd_key2data_key[[iName]]]][["penwidth"]] <- graph_nel@renderInfo@edges[["lwd"]][[iName]]
      #
      # graph_nel@renderInfo@edges[["arrowhead"]][[iName]] in theory should store the type of arrowhead
      # should be able to use the setter edgeData(graph, from, to, attr)
      # we need also to set the arrowtail
      # So maybe loading graph_nel and set arrowtail and arrowhead according to the amat_pag will be enough
      # this strength.plot code
      # https://github.com/cran/bnlearn/blob/414301e1a241148ec7bfb7e06f5eeda00ef2cd2b/R/frontend-plot.R#L34
      # calls inside graphviz.backend that manipulates arrowtail and arrowhead
      # https://github.com/cran/bnlearn/blob/414301e1a241148ec7bfb7e06f5eeda00ef2cd2b/R/graphviz-backend.R#L3
      # that inside plot with Rgraphviz::renderGraph(graph.plot)
      # maybe we just need to set arrowhead and tails and then a mod version of strength.plot calling a mod version of graphviz.backend that does not touch arrotail and arrowhead
      #
      # # cat(c(graph_nel@edgeData@data[[lwd_key2data_key[[iName]]]][["penwidth"]], "\n"))
    }
    # TODO fix it
    # toDot(graph = graph_nel, filename =  sub('.pdf', '.gv', tmp_fig))  # from library(Rgraphviz)
    #
    # TODO make it store edge weight
    #  toDotR(G = graph_nel, outDotFile = "~/Documents/figures/estonian.gv")  # from library(graph)
    # toDotR(G = graph_nel, outDotFile = sub('.pdf', '.gv', tmp_fig))
    # https://github.com/Bioconductor/graph/blob/master/R/TODOT.R
  }

  # FG Used when we were enforcing balacklist at the end
  # with_blacklist_avg_bn_obj <- avg_bn_obj
  # with_blacklist_bn_strength_obj <- bn_strength_obj
  # #
  # if ((original_bootstrap_n >= 1) & make_after_removal_plot_b) {
  #   cat(dir_arcs_to_be_removed_n, "directed edges must be removed", "\n")
  #   for (iPos_n in 1:dir_arcs_to_be_removed_n) {
  #     # Remove from graph bn object
  #     cat("Remove ",  dir_arcs_in_blacklist_mat[iPos_n, "from"], " -> ",  dir_arcs_in_blacklist_mat[iPos_n, "to"], "\n")
  #
  #     if (!show_graph_with_dashed_removed_b) {
  #       with_blacklist_avg_bn_obj <- drop.arc(x = with_blacklist_avg_bn_obj, from = dir_arcs_in_blacklist_mat[iPos_n, "from"], to = dir_arcs_in_blacklist_mat[iPos_n, "to"], debug = TRUE)
  #     }
  #     # Remove from bn.strength object
  #     # with_blacklist_bn_strength_obj <- as(
  #     #   filter(with_blacklist_bn_strength_obj, !(from == dir_arcs_in_blacklist_mat[iPos_n, "from"] & to == dir_arcs_in_blacklist_mat[iPos_n, "to"])),
  #     #   "bn.strength"
  #     # )
  #     #
  #     # with_blacklist_bn_strength_obj <- filter(with_blacklist_bn_strength_obj, !(from == dir_arcs_in_blacklist_mat[iPos_n, "from"] & to == dir_arcs_in_blacklist_mat[iPos_n, "to"]))
  #     # TODO vectorize it. Maybe with a group_by? In Pandas I would set index on ["from", "to"], create and index from dir_arcs_in_blacklist_mat
  #     with_blacklist_bn_strength_obj[
  #       (with_blacklist_bn_strength_obj$from == dir_arcs_in_blacklist_mat[iPos_n, "from"]) & (with_blacklist_bn_strength_obj$to == dir_arcs_in_blacklist_mat[iPos_n, "to"]),
  #       c("strength", "direction")
  #       ] <- 0
  #     #
  #     # drop.arc(x = bn_strength_obj, from = dir_arcs_in_blacklist_mat[iPos_n, "from"], to = dir_arcs_in_blacklist_mat[iPos_n, "to"], debug = TRUE)
  #   }
  #   attr(with_blacklist_bn_strength_obj, 'class') <- attr(bn_strength_obj, 'class')  # Because it lost the bn.stregth class attributes during manipulation
  #   # MEMO: if I use avg_bn_obj instead, the removed edges are present but dashed
  #   # if (show_graph_after_bgk_application_b) {
  #   plot_2 <- strength.plot(x = with_blacklist_avg_bn_obj, strength = with_blacklist_bn_strength_obj, threshold = attributes(bn_strength_obj)$threshold,
  #     main = second_title
  #   )
  #   # }
  # }
  # if ((original_bootstrap_n < 1) & make_after_removal_plot_b) {
  #   cat(dir_arcs_to_be_removed_n, "directed edges must be removed", "\n")
  #   for (iPos_n in 1:dir_arcs_to_be_removed_n) {
  #     # Remove from graph bn object
  #     cat("Remove ", dir_arcs_in_blacklist_mat[iPos_n, "from"], " -> ",  dir_arcs_in_blacklist_mat[iPos_n, "to"], "\n")
  #     bn_obj <- drop.arc(x = bn_obj, from = dir_arcs_in_blacklist_mat[iPos_n, "from"], to = dir_arcs_in_blacklist_mat[iPos_n, "to"], debug = TRUE)
  #   }
  #   plot_2 <- plot(bn_obj, main = second_title)
  # }

  # Failed attempt to use ggsave to save SVG
  # if (make_after_removal_plot_b) {
  #   grid.arrange(
  #     list(as.grob(plot_1), as.grob(plot_2)),
  #     nrow = show_graph_before_bgk_application_b + make_after_removal_plot_b
  #   )
  # } else {
  #   grid.arrange(as.grob(plot_1), nrow = show_graph_before_bgk_application_b + make_after_removal_plot_b)
  # }
  # #
  # path_to_fig <- config$get(option = "output_path_fig_pdf", section = "output_paths")
  # ggsave(
  #   filename = path_to_fig,
  #   width = 6, height = 3, onefile = FALSE,
  #   device = file_ext()
  # )
  dev.off()
  print_and_append_to_log(c("Completed", "\n"), fileConn)
  if (messaging_b) {
    sink(type = "message")  # reset message and output sink
    close(fileConn)
  }
}