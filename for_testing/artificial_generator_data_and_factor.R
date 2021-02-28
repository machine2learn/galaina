######################################################################################
# Generate random graph and thne data from those
# Author: Fabio Gori
######################################################################################


#### 0. Dependencies and Parameters ####

## dependent packages
library(pcalg)
library(infotheo)
library(BDgraph)
library(sbgcop)

library(comprehenr)
# source('for_testing/rmvCopulaFactorModel.R')

update_random_seed <- function(random_seed, parameter = 10) {
  random_seed + parameter
}

rmvCopulaFactorModel <- function(
  n = 1000, g = NULL, sigma = NULL,
  impurity = 0, range.indicators = 3:10, lambda.min = 0.1, lambda.max = 1, sd.residual = 1,
  first_random_seed = 1000,
  random_seed_update_parameter = 10
) {
  # This function aims to generate data from a Gaussian copula factor model
  #   given a graph or a covariance matrix in the latent space.
  #
  # Args:
  #   n, sample size
  #   g, a DAG (graphNEL)
  #   sigma, a covariance matrix
  #   impurity, No. of impurities of the measurement model
  #   range.indicators, range of the number of indicators per factor
  #   lambda.min and lambda.max, the minimum and maximum of factor loadings
  #   sd.residual, standard deviation of residuals
  #
  # Return: a list, which contains generated data and ...
  #

  library(mvtnorm)

  #### 0. Check inputs and Initialization
  ## check input g and sigma
  if (!is.null(g)) {
    ## pl: the number of latent factors (or variables)
    pl <- length(g@nodes)
    ## sigma: the population covariance matrix
    sigma <- trueCov(g)
  } else if (!is.null(sigma)) {
    pl <- ncol(sigma)
  } else {
    stop('Input g and sigma connot be empty simutanously.')
  }

  random_seed_n <- first_random_seed  # FG for controlling random sampling

  ## choose randomly some factors from pl as factors with more (>1) response variables
  # the number of factors with more (>1) response variables
  pm <- round(pl)
  index_more <- 1:pl  #sample(1:pl, pm)
  ## choose the number of response variables for each element of index_more from unif(3,10)
  #num_res_var = sample(3:10, pm, replace = T)#rep(4, pm)#
  if (length(range.indicators) == 1) {
    num_res_var <- rep(range.indicators, pm)
  } else {
    set.seed(random_seed_n)  # FG set random generator seed
    num_res_var <- sample(range.indicators, pm, replace = T)
    random_seed_n <- update_random_seed(random_seed_n, random_seed_update_parameter)  # FG update random seed
  }

  ## pr: the number of response variables.
  pr <- sum(num_res_var)


  #### 1. Generate normal data in latent sapce
  ## d_gauss stands for \eta
  set.seed(random_seed_n)  # FG set random generator seed
  d_gauss <- rmvnorm(n, sigma = sigma)#rmvDAG(n,g)
  random_seed_n <- update_random_seed(random_seed_n, random_seed_update_parameter)  # FG update random seed

  #### 2. Generate response data
  ## initializing
  Y <- mat.or.vec(n,pl+pr)
  Y[, 1:pl] <- d_gauss

  ## randomly generate factor loadings (Lamda) from latent factors to response variables
  Lamda <- matrix(0,pr,pm)
  index_res_var <- mat.or.vec(pm,1)
  for (i in 1:pm) {
    index_res_var[i] <- sum(num_res_var[1:i])
  }
  index_res_var <- c(0,index_res_var)
  for (i in 1:pm) {
    set.seed(random_seed_n)  # FG set random generator seed
    Lamda[(index_res_var[i]+1):index_res_var[i+1],i] <- runif(num_res_var[i], min = lambda.min, max = lambda.max)
    random_seed_n <- update_random_seed(random_seed_n, random_seed_update_parameter)  # FG update random seed
  }
  Lamda.pure <- Lamda
  ## noise
  set.seed(random_seed_n)  # FG set random generator seed
  error <- matrix(rnorm(pr*n, 0, sd.residual),n,pr)
  random_seed_n <- update_random_seed(random_seed_n, random_seed_update_parameter)  # FG update random seed
  ## add some impurities
  if (impurity > 0){
    # the number of type 1 and 2 impurities
    ni.1 <- round(impurity)
    ni.2 <- impurity - ni.1
    # type 1: to Lamda
    if (ni.1 > 0) {
      set.seed(random_seed_n)  # FG set random generator seed
      Lamda[sample(which(Lamda == 0), ni.1)] <- runif(ni.1, 0.1, 1)
      random_seed_n <- update_random_seed(random_seed_n, random_seed_update_parameter)  # FG update random seed
    }
    # type 2: to error
    if (ni.2 > 0){
      vec <- rep(0, pr*(pr-1)/2)
      set.seed(random_seed_n)  # FG set random generator seed
      vec[sample(pr*(pr-1)/2, ni.2)] <- runif(ni.2)
      random_seed_n <- update_random_seed(random_seed_n, random_seed_update_parameter)  # FG update random seed
      sigma.error <- diag(pr)
      sigma.error[upper.tri(sigma.error)] <- vec
      sigma.error[lower.tri(sigma.error)] <- t(sigma.error)[lower.tri(sigma.error)]
      set.seed(random_seed_n)  # FG set random generator seed
      error <- rmvnorm(n, sigma = sigma.error)
      random_seed_n <- update_random_seed(random_seed_n, random_seed_update_parameter)  # FG update random seed
    }
  }


  ## generate responsing data
  Y[,(pl+1):(pl+pr)] <- d_gauss[,index_more] %*% t(Lamda) + error


  #### 3. Return
  list(data = Y[,-(1:pl)], Sigma = sigma, index_more = index_more, Lambda = Lamda.pure)
}

single_artificial_generator_data_and_factor <- function(
  factor_n = 7,
  prob = NULL,
  lB = 0.1, uB = 1,
  frac_variable_discreet = 0.5,
  range_bin_discreet = 2:5,
  #rmv parameters
  n = 1000,
  g = NULL,
  impurity = 0, range.indicators = 3:10, lambda.min = 0.1, lambda.max = 1, sd.residual = 1,
  first_random_seed = 1000,
  random_seed_update_parameter = 10
) {
  ## parameters
  # No. of factors
  # factor_n <- 6
  # sample size
  # n <- 1000
  if (is.null(prob)) {
    prob <- 2 / (factor_n - 1)
  }

  #### 1. Simulate Data ####

  ## simulate a DAG
  # a random DAG, serving as the true DAG in latent space
  random_seed_n <- first_random_seed  # FG for controlling random sampling
  set.seed(random_seed_n)  # FG set random generator seed
  g <- randomDAG(n = factor_n, prob = prob, lB = lB, uB = uB)
  random_seed_n <- update_random_seed(random_seed_n, random_seed_update_parameter)  # FG update random seed

  # true CPDAG
  # g.cpdag <- dag2cpdag(g)
  ## generate data
  set.seed(random_seed_n)  # FG set random generator seed

  data.obs <- rmvCopulaFactorModel(
    n = n,
    g = g,
    impurity = impurity,
    range.indicators = range.indicators,
    lambda.min = lambda.min,
    lambda.max = lambda.max,
    sd.residual = sd.residual,
    first_random_seed = random_seed_n,
    random_seed_update_parameter = random_seed_update_parameter
  )
  random_seed_n <- update_random_seed(random_seed_n, random_seed_update_parameter)  # FG update random seed
  # mapping from latents to response variables
  Lambda <- (data.obs$Lambda != 0) * 1
  # data in response space (fully Gaussian)
  Z <- data.obs$data
  # data in observed space (mixed continuous and ordinal)
  var_n <- ncol(Z)
  Y <- Z
  set.seed(random_seed_n)  # FG set random generator seed
  sampleResult <- sample(x = 1:var_n, size = round(var_n * frac_variable_discreet))
  random_seed_n <- update_random_seed(random_seed_n, random_seed_update_parameter)  # FG update random seed
  for (i in sampleResult) {
    set.seed(random_seed_n)  # FG set random generator seed
    Y[, i] <- matrix(unlist(discretize(Z[, i], nbins = sample(range_bin_discreet, 1))), byrow = FALSE, nrow = n)
    random_seed_n <- update_random_seed(random_seed_n, random_seed_update_parameter)  # FG update random seed
  }
  out_ls <- list(
    original_graph = g,
    for_inference = list(Y = Y, Lambda = Lambda),
    generator_parameters = list(
      factor_n = factor_n,
      prob = prob,
      lB = lB,
      uB = uB,
      frac_variable_discreet = frac_variable_discreet,
      range_bin_discreet = range_bin_discreet,
      n = n,
      g = g,
      impurity = impurity,
      range.indicators = range.indicators,
      lambda.min = lambda.min,
      lambda.max = lambda.max,
      sd.residual = sd.residual,
      first_random_seed = first_random_seed,
      random_seed_update_parameter = random_seed_update_parameter
    )
  )
  return(out_ls)

  # #### 2. Inference ####
  #
  # ## inference
  # cop.fac.obj <- inferCopulaFactorModel(Y, Lambda = Lambda, nsamp = 1000)
  # # extract samples of the correlation matrix over latent variables, ignoring the first 500 samples (burn-in)
  # C.samples <- cop.fac.obj$Sigma.psamp[1:factor_n, 1:factor_n, 501:1000]
  # # the posterior mean
  # C <- apply(C.samples, c(1, 2), mean)
  # # standard deviations
  # C.sd <- apply(C.samples, c(1, 2), sd)
  # # effective sample size
  # C.ess <- ((1 - C^2) / C.sd)^2
  #
  # #### 3. Causal discovery ####
  #
  # ## call the order independent version of the standard PC algorithm
  # graph.cfpc <- pc(suffStat = list(C = C, n = mean(C.ess[upper.tri(C.ess)])), indepTest = gaussCItest,
  #                  alpha = 0.05, p = factor_n, skel.method = "stable", maj.rule = T, solve.confl = T)
  #
  # ## show results
  # par(mfrow = c(1, 2))
  # plot(g.cpdag, main = 'True Graph')
  # plot(graph.cfpc, main = 'Copula Factor PC')
}

multi_artificial_generator_data_and_factor <- function(
  dataset_n = 1,
  factor_n = 7,
  prob = NULL,
  lB = 0.1, uB = 1,
  frac_variable_discreet = 0.5,
  range_bin_discreet = 2:5,
  #rmv parameters
  n = 1000,
  g = NULL,
  impurity = 0, range.indicators = 3:10, lambda.min = 0.1, lambda.max = 1, sd.residual = 1,
  first_random_seed = 1000,
  random_seed_update_parameter = 10
) {
  multi_out_ls <- list()
  random_seed_n <- first_random_seed  # FG for controlling random sampling
  for (i_run in 1:dataset_n) {
    set.seed(random_seed_n)
    tmp_ls <- single_artificial_generator_data_and_factor(
      factor_n = factor_n,
      prob = prob,
      lB = lB,
      uB = uB,
      frac_variable_discreet = frac_variable_discreet,
      range_bin_discreet = range_bin_discreet,
      n = n,
      g = g,
      impurity = impurity,
      range.indicators = range.indicators,
      lambda.min = lambda.min,
      lambda.max = lambda.max,
      sd.residual = sd.residual,
      first_random_seed = random_seed_n,
      random_seed_update_parameter = random_seed_update_parameter
    )
    random_seed_n <- update_random_seed(random_seed_n, random_seed_update_parameter)  # FG update random seed
    multi_out_ls[[i_run]] <- tmp_ls
  }
  return(multi_out_ls)
}