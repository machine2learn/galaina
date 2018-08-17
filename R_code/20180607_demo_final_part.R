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

#
factor_n <- length(factor_names)

#### 2. Perform inference ####
## some checks
stopifnot(all(sapply(data_df, is.numeric)))  # check all columns in df are numeric
stopifnot(all(sapply(factor_loading_df, is.numeric)))  # check all columns in df are numeric

cop_fac_obj <- my_inferCopulaFactorModel(Y = data.matrix(data_df), Lambda = data.matrix(factor_loading_df), nsamp = gibs_sampling_n)  # , tol = 1e-22)
# cop_fac_obj <- my_inferCopulaFactorModel(Y = as.matrix(data_df), Lambda = factor_loading_df, nsamp = gibs_sampling_n, tol = 1e-22)
# cop_fac_obj <- inferCopulaFactorModel(Y = data_matrix, Lambda = factor_loading_matrix, nsamp = gibs_sampling_n)
# extract samples of the correlation matrix over latent variables, ignoring the first 500 samples (burn-in)
C_samples <- cop_fac_obj$Sigma.psamp[1:factor_n, 1:factor_n, (throw_away_n + 1) : gibs_sampling_n]
# the posterior mean
C <- apply(C_samples, c(1,2), mean)

C_sd <- apply(C_samples, c(1,2), sd)
# effective sample size
C_ess <- ((1-C^2)/C_sd)^2

## export to CSV
write.csv(C, file = args[6])

#### 3. Causal discovery ####

## call the order independent version of the standard PC algorithm
graph_cfpc <- pc(
  suffStat = list(C = C, n = mean(C_ess[upper.tri(C_ess)])), indepTest = gaussCItest, 
  labels = factor_names,
  alpha = 0.05, skel.method = "stable", maj.rule = TRUE, solve.confl = TRUE
)
saveRDS(graph_cfpc, file = args[7])

## show results
graphics.off()
# pdf('demo_aggressotype.pdf', width=6, height=3)
pdf(args[3], width=6, height=3)
# plot.new()
# par(mfrow = c(1,2))
# par(mfrow = c(1, 1))
# plot(g.cpdag, main = 'True Graph')
plot(graph_cfpc, main = 'Copula Factor PC')
dev.off()
# print(length(factor_names))

# graph_cfpc <- readRDS(file = "test_consistent_55_PcAlgo.rds")
cfpc_bn <- as.bn(x = graph_cfpc)