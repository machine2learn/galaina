#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(utils)
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
data_df <- read.csv(file = args[1], row.names = 1, na.strings=c("",".","NA"))
factor_loading_df <- read.csv(file = args[2], row.names = 1,  na.strings=c("",".","NA"))

factor_names <- colnames(factor_loading_df)
# load(file = '/Volumes/encrypted_preprocessed_data/demo_aggressotype_transfer.rda')
# data_matrix <- as.numeric(data_matrix)
# factor_loading_matrix <- as.numeric(factor_loading_matrix)

# print(typeof(data_matrix))
# print(typeof(factor_loading_matrix))

gibs_sampling_n <- as.integer(args[4])  # 1000

throw_away_n <- as.integer(args[5])  # 500
bootstrap_n <-  as.integer(args[6])  # 500
cat("Gibbs sampling: ", gibs_sampling_n, "\n")
cat("Throw away: ", throw_away_n, "\n")
cat("Bootstrap: ", bootstrap_n, "\n")

first_random_seed_n <- 100  # TODO make it a parameter
random_seed_update_parameter <- 10 # TODO see if it should be a parameter
causal_discovery_observation_n <- 0  # TODO make it a parameter

odens_n <- 1  # save only every odens_n 
#
factor_n <- length(factor_names)

original_first_random_seed_n <- first_random_seed_n
random_seed_n <- first_random_seed_n

#### 2. Perform bootstrap, inference, causal discovery ####
## some checks
stopifnot(throw_away_n < gibs_sampling_n)
stopifnot(all(sapply(data_df, is.numeric)))  # check all columns in df are numeric
stopifnot(all(sapply(factor_loading_df, is.numeric)))  # check all columns in df are numeric

bootstrapped_pcalgo_list <- vector(mode = "list", length = bootstrap_n)
throw_away_odens_n <- floor(throw_away_n / odens_n)
for (iBootstrap_n in 1:bootstrap_n) {
  ## 2.1 Bootstrap
  cat("Bootstrap run: ", iBootstrap_n, "\n")
  set.seed(random_seed_n)
  tmp_bootstrap_data_df <- data_df[sample(nrow(data_df)), ]
  random_seed_n <- update_random_seed(random_seed_n, random_seed_update_parameter)  # FG update random seed
  
  ## 2.2 Perform inference
  tmp_cop_fac_obj <- my_inferCopulaFactorModel(
    Y = data.matrix(tmp_bootstrap_data_df), 
    Lambda = data.matrix(factor_loading_df), 
    nsamp = gibs_sampling_n,
    odens = odens_n,  # Store each Gibbs sampling output
    first_random_seed = original_first_random_seed_n,
    random_seed_update_parameter = random_seed_update_parameter
  )  # , tol = 1e-22)
  
  cat(dim(tmp_cop_fac_obj$Sigma.psamp), "\n")
  # extract samples of the correlation matrix over latent variables, ignoring the first samples (burn-in)
  # C_samples <- tmp_cop_fac_obj$Sigma.psamp[1:factor_n, 1:factor_n, (throw_away_n + 1) : gibs_sampling_n]
  if (throw_away_odens_n < 1) {
    C_samples <- tmp_cop_fac_obj$Sigma.psamp[1:factor_n, 1:factor_n, ]
  } else {
    C_samples <- tmp_cop_fac_obj$Sigma.psamp[1:factor_n, 1:factor_n, -(1:throw_away_odens_n)]
  }
  # the posterior mean
  C <- apply(C_samples, c(1,2), mean)
  if (causal_discovery_observation_n == 0) {
    # TODO change to allow user to input a fix value to C_ess
    C_sd <- apply(C_samples, c(1,2), sd)
    # effective sample size
    C_ess <- ((1-C^2)/C_sd)^2
    # effective_observation_n <- mean(C_ess[upper.tri(C_ess)])
    causal_discovery_observation_n <- mean(C_ess[upper.tri(C_ess)])
  }
  
  ## export to CSV
  # write.csv(C, file = args[7])
  
  #### 2.3 Causal discovery ####
  # MEMO we could directly store the bnlearn version
  ## call the order independent version of the standard PC algorithm
  # tmp_graph_cfpc <- pc(
  bootstrapped_pcalgo_list[[iBootstrap_n]] <- pc(
    suffStat = list(
      C = C,
      n = causal_discovery_observation_n  # mean(C_ess[upper.tri(C_ess)])
    ), 
    indepTest = gaussCItest, 
    labels = factor_names,
    alpha = 0.05, 
    skel.method = "stable", 
    maj.rule = TRUE, 
    solve.confl = TRUE
  )
  # bootstrapped_graph_list[[iBootstrap_n]] <- tmp_graph_cfpc
  # saveRDS(graph_cfpc, file = args[8])
}

saveRDS(bootstrapped_pcalgo_list, file = args[8])
bootstrapped_bnlearn_list <- lapply(X = bootstrapped_pcalgo_list, FUN = as.bn)
bn_strength_obj <- custom.strength(networks = bootstrapped_bnlearn_list, nodes = factor_names)
saveRDS(bn_strength_obj, file = "/Users/fabio/bn_strength_obj.rds")

avg_bn_obj <- averaged.network(strength = bn_strength_obj, nodes = factor_names)  # it automatically gets the thresold from bn_strength_obj
saveRDS(avg_bn_obj, file = "/Users/fabio/avg_bn_obj.rds")

graphics.off()
pdf(args[3], width=6, height=3)
strength.plot(x = avg_bn_obj, strength = bn_strength_obj, threshold = attributes(bn_strength_obj)$threshold)
dev.off()

# graph_cfpc <- readRDS(file = "test_consistent_55_PcAlgo.rds")
# cfpc_bn <- as.bn(x = graph_cfpc)

# ## show results
# graphics.off()
# pdf(args[3], width=6, height=3)
# plot(graph_cfpc, main = 'Copula Factor PC')
# dev.off()