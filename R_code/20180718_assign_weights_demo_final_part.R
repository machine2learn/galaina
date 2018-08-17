#!/usr/bin/env Rscript

# MAIN

# args = commandArgs(trailingOnly=TRUE)
args <- c("/Volumes/encrypted_data_sparse/demo_output/test_background_1000_bootstrap_3_merged_data.csv",
          "/Volumes/encrypted_data_sparse/demo_output/test_background_1000_bootstrap_3_merged_factor_loading_matrix.csv",
          "/Volumes/encrypted_data_sparse/demo_output/test_background_1000_bootstrap_1_causal_discovery.pdf",
          "1000", "500", "0",
          # "/Volumes/encrypted_data_sparse/demo_output/test_background_1000_bootstrap_3_correlation_matrix.csv",
          "/Volumes/encrypted_data_sparse/demo_output/test_background_1000_bootstrap_1_PcAlgo.rds",
          "/Volumes/encrypted_data_sparse/data_for_demo/background_knoledge.csv",
          "0"
          # "70"  # we need to store here the sample size
          )
          # "/Volumes/encrypted_data_sparse/data_for_demo/background_knoledge_exp.csv",  # edited edge
column_separator_str = ',' <- "NULL"
library(utils)
library(dplyr)  # actually I load it only for filter
library(stringi)  # to check if strings are empty
library(sbgcop)
library(infotheo)
library(BDgraph)
library(pcalg)  # depends on bioconductor
library(bnlearn)

# source('/Users/fabio/projects/aggressotype/external_code/CopulaFactorModel/R/my_inferCopulaFactorModel.R')
source('/Users/fabio/projects/aggressotype/code/software/software_demo/R_code/my_inferCopulaFactorModel.R')
# source('my_inferCopulaFactorModel.R')  # TODO test it

#### 1. Load input ####
# Load data_matrix, factor_loading_matrix
# data_df <- read.csv(file = '/Volumes/encrypted_preprocessed_data/data_for_demo/intersection/data_df.csv', row.names = 1, na.strings=c("",".","NA"))
# factor_loading_df <- read.csv(file = '/Volumes/encrypted_preprocessed_data/data_for_demo/intersection/factor_loading_df.csv', row.names = 1,  na.strings=c("",".","NA"))


data_df <- read.csv(file = args[1], sep = column_separator_str, row.names = 1, na.strings = c("", ".", "NA"))

factor_loading_df <- read.csv(file = args[2], sep = column_separator_str, row.names = 1, na.strings = c("", ".", "NA"))

factor_names <- colnames(factor_loading_df)
# load(file = '/Volumes/encrypted_preprocessed_data/demo_aggressotype_transfer.rda')
# data_matrix <- as.numeric(data_matrix)
# factor_loading_matrix <- as.numeric(factor_loading_matrix)

# print(typeof(data_matrix))
# print(typeof(factor_loading_matrix))

gibs_sampling_n <- as.integer(args[4])  # in INI e.g. 1000 
gibbs_burn_in_n <- as.integer(args[5])  # in INI e.g. 500
#
gibbs_first_random_seed_n <- 365  # TODO make it a option in INI
gibbs_random_seed_update_parameter_n <- 10  # TODO see if it should be a option in INI

original_bootstrap_n <- as.integer(args[6])  # in INI e.g. 500
bootstrap_n < max(original_bootstrap_n, 1)
#
bootstrap_first_random_seed_n <- 100  # TODO make it a option in INI
bootstrap_random_seed_update_parameter_n <- 10  # TODO see if it should be a option in INI

input_path_directed_edges_blacklist <- args[8]

cat("Gibbs total samples:", gibs_sampling_n, "\n")
cat("Gibbs burn-in samples:", gibbs_burn_in_n, "\n")
# cat("Bootstrap samples: ", bootstrap_n, "\n")

# TODO add input matrix with two columns: "from", to" - or "x", "y"
# y <- c('COR_1', 'COR_10', 'COR_11', 'COR_12', 'COR_13', 'COR_14', 'COR_15',
#   'COR_16', 'COR_17', 'COR_18', 'COR_19', 'COR_2', 'COR_20', 'COR_21',
#   'COR_3', 'COR_4', 'COR_5', 'COR_6', 'COR_7', 'COR_8', 'COR_9'
# )
# x <- rep(c('Gender'), times = length(y))
# #
# input_bgk_from_to_mat <- cbind(from=x, to=y)

# TODO future: use whitelist and blacklist
if (file.exists(input_path_directed_edges_blacklist)) {
  dir_arc_blacklist_mat <- as.matrix(read.csv(file = input_path_directed_edges_blacklist, sep = column_separator_str, row.names = NULL))
  # input_bgk_whitelist_from_to_mat <- as.matrix(read.csv(file = args[9], sep = column_separator_str, row.names = NULL))
}
# if (file.exists(args[10])) {
#   input_bgk_blacklist_from_to_mat <- as.matrix(read.csv(file = args[10], sep = column_separator_str, row.names = NULL))
# }
# input_bgk_from_to_mat <- cbind(y, x)  # better use form above because this is to from and now we useded the ones with rev columns at the same time
##

causal_discovery_observation_n <- as.integer(args[9])  # TODO make it a option in INI

odens_n <- 1  # save only every odens_n
#
factor_n <- length(factor_names)

bootstrap_random_seed_n <- bootstrap_first_random_seed_n

#### 2. Perform bootstrap, inference, causal discovery ####
## some checks
stopifnot(gibbs_burn_in_n < gibs_sampling_n)
stopifnot(all(sapply(data_df, is.numeric)))  # check all columns in df are numeric
stopifnot(all(sapply(factor_loading_df, is.numeric)))  # check all columns in df are numeric

bootstrapped_pcalgo_list <- vector(mode = "list", length = bootstrap_n)
# bootstrapped_graphnel_list <- vector(mode = "list", length = bootstrap_n)
throw_away_odens_n <- floor(gibbs_burn_in_n / odens_n)

for (iBootstrap_n in 1:bootstrap_n) {
  ## 2.1 Bootstrap
  if (original_bootstrap_n == 0) {
    cat("No bootstrap performed \n")
    tmp_bootstrap_data_df <- data_df
  } else {
    cat("Current bootstrap sample:", iBootstrap_n, "of", bootstrap_n, "\n")
    set.seed(bootstrap_random_seed_n)
    tmp_bootstrap_data_df <- data_df[sample(nrow(data_df)), ]
    bootstrap_random_seed_n <- update_random_seed(bootstrap_random_seed_n, bootstrap_random_seed_update_parameter_n)  # FG update random seed
  }
 
  ## 2.2 Perform inference
  tmp_cop_fac_obj <- my_inferCopulaFactorModel(
    Y = data.matrix(tmp_bootstrap_data_df),
    Lambda = data.matrix(factor_loading_df),
    nsamp = gibs_sampling_n,
    odens = odens_n,  # Store each Gibbs sampling output
    first_random_seed = gibbs_first_random_seed_n,
    random_seed_update_parameter = gibbs_random_seed_update_parameter_n
  )  # , tol = 1e-22)
 
  # cat(dim(tmp_cop_fac_obj$Sigma.psamp), "\n")
  # extract samples of the correlation matrix over latent variables, ignoring the first samples (burn-in)
  # C_samples <- tmp_cop_fac_obj$Sigma.psamp[1:factor_n, 1:factor_n, (gibbs_burn_in_n + 1) : gibs_sampling_n]
  if (throw_away_odens_n < 1) {
    C_samples <- tmp_cop_fac_obj$Sigma.psamp[1:factor_n, 1:factor_n, ]
  } else {
    C_samples <- tmp_cop_fac_obj$Sigma.psamp[1:factor_n, 1:factor_n, -(1:throw_away_odens_n)]
  }
  # the posterior mean
  C <- apply(C_samples, c(1,2), mean)
  if (causal_discovery_observation_n == 0) {
    C_sd <- apply(C_samples, c(1,2), sd)
    # TODO change to allow user to input a fix value to causal_discovery_observation_n
    # effective sample size
    C_ess <- ( (1 - C^2) / C_sd )^2
    # effective_observation_n <- mean(C_ess[upper.tri(C_ess)])
    causal_discovery_observation_n <- mean(C_ess[upper.tri(C_ess)])
  }

  ## export to CSV
  # write.csv(C, file = args[7])

  # set.seed(random_seed_n)
  #
  #### 2.3 Causal discovery ####
  # MEMO we could directly store the bnlearn version
  ## call the order independent version of the standard PC algorithm
  # bootstrapped_pcalgo_list[[iBootstrap_n]] <- pc(
  tmp_graph_cfpc <- pc(
    suffStat = list(
      C = C,
      n = causal_discovery_observation_n  # mean(C_ess[upper.tri(C_ess)])
    ),
    indepTest = gaussCItest,
    alpha = 0.05,
    labels = factor_names,
    fixedGaps = NULL,
    fixedEdges = NULL,
    NAdelete = TRUE,
    m.max = Inf,
    u2pd = "relaxed",  # impled by solve.confl = TRUE
    skel.method = "stable",
    # conservative = TRUE,
    # Choose either conservative PC or majority rule PC!
    maj.rule = TRUE,
    solve.confl = TRUE
  )
  # random_seed_n <- update_random_seed(random_seed_n, random_seed_update_parameter_n)  # FG update random seed
  #
  # WARNING we might have to turn off the edge assignment part
 
  # tmp_amat_cpdag <- as(tmp_graph_cfpc, "amat")
  # with_bgk_amat_cpdag <- addBgKnowledge(gInput = tmp_amat_cpdag, checkInput = FALSE)  # apply orientation rules of Meek 1995
  # DEBUG begin
  # cat(is.matrix(tmp_amat_cpdag), "\n")

  # DEBUG begin
  # Restrict to arc that can be enforced on this graph: the must correspond to egdes that are undirected
  # Check if edge is there
 
  # if (file.exists(args[9])) {
  #   tmp_arc_enforceable_bv <- as.logical(tmp_amat_cpdag[input_bgk_from_to_mat]) & as.logical(tmp_amat_cpdag[input_bgk_from_to_mat[, rev(colnames(input_bgk_from_to_mat))]])
  #   # tmp_arc_enforceable_bv <- as.logical(tmp_amat_cpdag[input_bgk_from_to_mat])
  #   # cat(tmp_arc_enforceable_bv, "\n")
  #   # If some enforceable arc are still there, enforce them
  #   if (any(tmp_arc_enforceable_bv)) {
  #     # Maybe check x and y are not empty
  #     with_bgk_amat_cpdag <- addBgKnowledge(
  #       gInput = tmp_amat_cpdag,
  #       # x = x[tmp_arc_enforceable_bv],
  #       # y = y[tmp_arc_enforceable_bv],
  #       x = input_bgk_from_to_mat[tmp_arc_enforceable_bv, "from"],
  #       y = input_bgk_from_to_mat[tmp_arc_enforceable_bv, "to"],
  #       verbose = TRUE,
  #       checkInput = FALSE  # because it was giving error
  #     )
  #     # cat(with_bgk_amat_cpdag, "\n")
  #     # cat(typeof(with_bgk_amat_cpdag), "\n")
  #   }
  # }
  # bootstrapped_graphnel_list[[iBootstrap_n]] <- as(t(as(with_bgk_amat_cpdag, "matrix")), "graphNEL")  # strange we need to specify matrix
 
  # In theory we would need to convert to grphNEL only if we add background knowledge, because as.bn would work with both
  # bootstrapped_graphnel_list[[iBootstrap_n]] <- as(t(with_bgk_amat_cpdag), "graphNEL")
  bootstrapped_pcalgo_list[[iBootstrap_n]] <- tmp_graph_cfpc
  # saveRDS(graph_cfpc, file = args[8])
  cat('\n\n')
}
cat('\n')

# TODO adjust it when no bootstap is done
if (!stri_isempty(args[8])) {
  saveRDS(bootstrapped_pcalgo_list, file = args[8])
}
if (original_bootstrap_n > 0) {
  bootstrapped_bnlearn_list <- lapply(X = bootstrapped_pcalgo_list, FUN = as.bn)
  # saveRDS(bootstrapped_graphnel_list, file = args[8])
  # bootstrapped_bnlearn_list <- lapply(X = bootstrapped_graphnel_list, FUN = as.bn)
  bn_strength_obj <- custom.strength(networks = bootstrapped_bnlearn_list, nodes = factor_names)
  
  path_to_bn_strength_obj = "/Users/fabio/bn_strength_obj.rds"  # TODO put in INI
  if (!stri_isempty(path_to_bn_strength_obj)) {
    saveRDS(bn_strength_obj, file = path_to_bn_strength_obj)
  }
  
  avg_bn_obj <- averaged.network(strength = bn_strength_obj, nodes = factor_names)  # it automatically gets the thresold from bn_strength_obj
  # saveRDS(avg_bn_obj, file = "/Users/fabio/avg_bn_obj.rds")
  path_to_avg_bn_obj_before_bgk_application = "/Users/fabio/avg_bn_obj.rds"  # TODO put in INI
  if (!stri_isempty(path_to_avg_bn_obj_before_bgk_application)) {
    saveRDS(avg_bn_obj, file = path_to_avg_bn_obj_before_bgk_application)
  }
} else {
  bn_obj <- as.bn(bootstrapped_pcalgo_list[[1]])
}
  
main_plot_title_str <- "Fake model"  # TODO put in INI
show_graph_before_bgk_application_b <- as.logical("False")  # TODO put in INI
show_graph_after_bgk_application_b <- as.logical("True")  # TODO put in INI
show_graph_with_dashed_removed_b <- as.logical("False")  # TODO put in INI

# Check if there are some blacklist directed edges to be removed
if (original_bootstrap_n > 0){
  dir_arcs_mat <- directed.arcs(avg_bn_obj)
} else {
  dir_arcs_mat <- directed.arcs(bn_obj)
}
dir_arcs_to_be_removed_n <- 0

if ((nrows(dir_arc_mat) > 0) & file.exists(input_path_directed_edges_blacklist) & (nrows(dir_arc_blacklist_mat) > 0)) {  # if at least one od those mat is empty, there are no arcs to be removed
  dir_arcs_and_blacklist_mat <- rbind(unique(dir_arcs_mat), unique(dir_arc_blacklist_mat))
  dir_arcs_in_blacklist_mat <- dir_arcs_and_blacklist_mat[duplicated(dir_arcs_and_blacklist_mat), , drop = FALSE]
  dir_arcs_to_be_removed_n <- nrow(dir_arcs_in_blacklist_mat)
}

make_after_removal_plot_b <- (show_graph_after_bgk_application_b & (dir_arcs_to_be_removed_n > 0))
if (!(make_after_removal_plot_b | show_graph_before_bgk_application_b)) {
  show_graph_before_bgk_application_b <- TRUE  # if there is no second plot to make, force to make the first one
}

graphics.off()

pdf(args[3], width=6, height=3, onefile = FALSE)  # TODO take it from INI
par(mfrow=c(show_graph_before_bgk_application_b + make_after_removal_plot_b, 1))
par(mar=c(1,1,1,1))
first_title <- main_plot_title_str  # 'Causal graph with all edges'
if (make_after_removal_plot_b) {
  first_title <- paste(main_plot_title_str, "- including blacklist")  #'Before removing blacklist edges'
}
if (!(make_after_removal_plot_b) & (dir_arcs_to_be_removed_n == 0)) {
  first_title <- paste(main_plot_title_str, "- no blacklist effect")  # No edges were removed
}
second_title <- main_plot_title_str
if (make_after_removal_plot_b & show_graph_with_dashed_removed_b) {
  second_title <- paste(main_plot_title_str, "- blacklist dashed")  #'Before removing blacklist edges'
} else if (make_after_removal_plot_b & !show_graph_with_dashed_removed_b) {
  second_title <- paste(main_plot_title_str, "- excluding blacklist")  #'Before removing blacklist edges'
}

cat(first_title, "\n")
if (show_graph_before_bgk_application_b) {
  if (original_bootstrap_n > 0) {
    strength.plot(x = avg_bn_obj, strength = bn_strength_obj, threshold = attributes(bn_strength_obj)$threshold, main = first_title)
  } else {
    plot(bootstrapped_pcalgo_list[[1]], main = first_title)  # if 
  }
}
# TODO
# D < - get directed edges
# if (original_bootstrap_n > 0) {
with_blacklist_avg_bn_obj <- avg_bn_obj
with_blacklist_bn_strength_obj <- bn_strength_obj
#
if ((original_bootstrap_n > 0) & make_after_removal_plot_b) {
  cat(dir_arcs_to_be_removed_n, "directed edges must be removed", "\n")
  for (iPos_n in 1:dir_arcs_to_be_removed_n) {
    # Remove from graph bn object
    cat("Remove ",  dir_arcs_in_blacklist_mat[iPos_n, "from"], " -> ",  dir_arcs_in_blacklist_mat[iPos_n, "to"], "\n")
    
    if (!show_graph_with_dashed_removed_b) {
      with_blacklist_avg_bn_obj <- drop.arc(x = with_blacklist_avg_bn_obj, from = dir_arcs_in_blacklist_mat[iPos_n, "from"], to = dir_arcs_in_blacklist_mat[iPos_n, "to"], debug = TRUE)
    }
    # Remove from bn.strength object
    # with_blacklist_bn_strength_obj <- as(
    #   filter(with_blacklist_bn_strength_obj, !(from == dir_arcs_in_blacklist_mat[iPos_n, "from"] & to == dir_arcs_in_blacklist_mat[iPos_n, "to"])),
    #   "bn.strength"
    # )
    #
    # with_blacklist_bn_strength_obj <- filter(with_blacklist_bn_strength_obj, !(from == dir_arcs_in_blacklist_mat[iPos_n, "from"] & to == dir_arcs_in_blacklist_mat[iPos_n, "to"]))
    # TODO vectorize it. Maybe with a group_by? In Pandas I would set index on ["from", "to"], create and index from dir_arcs_in_blacklist_mat
    with_blacklist_bn_strength_obj[
      (with_blacklist_bn_strength_obj$from == dir_arcs_in_blacklist_mat[iPos_n, "from"]) & (with_blacklist_bn_strength_obj$to == dir_arcs_in_blacklist_mat[iPos_n, "to"]),
      c("strength", "direction")
      ] <- 0
    #
    # drop.arc(x = bn_strength_obj, from = dir_arcs_in_blacklist_mat[iPos_n, "from"], to = dir_arcs_in_blacklist_mat[iPos_n, "to"], debug = TRUE)
  }
  attr(with_blacklist_bn_strength_obj, 'class') <- attr(bn_strength_obj, 'class')  # Because it lost the bn.stregth class attributes during manipulation
  # MEMO: if I use avg_bn_obj instead, the removed edges are present but dashed
  # if (show_graph_after_bgk_application_b) {
  strength.plot(x = with_blacklist_avg_bn_obj, strength = with_blacklist_bn_strength_obj, threshold = attributes(bn_strength_obj)$threshold,
    main = second_title
  )
  # }
}
if ((original_bootstrap_n = 0) & make_after_removal_plot_b) {
  cat(dir_arcs_to_be_removed_n, "directed edges must be removed", "\n")
  for (iPos_n in 1:dir_arcs_to_be_removed_n) {
    # Remove from graph bn object
    cat("Remove ",  dir_arcs_in_blacklist_mat[iPos_n, "from"], " -> ",  dir_arcs_in_blacklist_mat[iPos_n, "to"], "\n")
    bn_obj <- drop.arc(x = bn_obj, from = dir_arcs_in_blacklist_mat[iPos_n, "from"], to = dir_arcs_in_blacklist_mat[iPos_n, "to"], debug = TRUE)
  }
  plot(bn_obj, main = second_title)
}
dev.off()
# -
# for i
# get undirected edges


# pdf(args[3], width=6, height=3)

# graph_cfpc <- readRDS(file = "test_consistent_55_PcAlgo.rds")
# cfpc_bn <- as.bn(x = graph_cfpc)

# ## show results
# graphics.off()
# pdf(args[3], width=6, height=3)
# plot(graph_cfpc, main = 'Copula Factor PC')
# dev.off()