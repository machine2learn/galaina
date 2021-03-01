# library(here)
# source(here("R_code/20180725_use_config_ini_final_part.R"))  # Load print_and_append_to_log
print_and_append_to_log <- function(tmp_log_str_v, fileConn = NULL) {
  cat(tmp_log_str_v)
  if (!is.null(fileConn) & !is.list(fileConn)) {  # I put this not-list confition because in this way it skips also fileConn = list(NULL)
    write(paste(tmp_log_str_v, collapse = " "), file = fileConn, append = TRUE)
  }
}

update_random_seed <- function(random_seed, parameter = 10) {
  random_seed + parameter
}

version_at_least <- function(pkg, than) {
  as.logical((compareVersion(as.character(packageVersion(pkg)), than) >= 0))
}

stopGeneric <- function(condition_b, fileConn, msg) {
  if (!(condition_b)) {
    # if (!is.null(fileConn)) {
    print_and_append_to_log(c(msg, "\n"), fileConn)
    # }
    stop(msg)
  }
}

showData <- function(focus_v, fileConn, msg) {
  print_and_append_to_log(c(msg, "\n", sep = ""), fileConn)
  print_and_append_to_log(c(capture.output(focus_v), "\n", sep = "\n"), fileConn)
}
# source("shared_functions.R")
# FG Keep these two here
# Y_is_finite <- all(is.finite(Y))
# X_is_finite <- all(is.finite(X))
this_should_happen <- function(Y, X) {
  Y_is_finite <- all(is.finite(Y))
  X_is_finite <- all(is.finite(X))
  as.logical(
    (Y_is_finite & X_is_finite) | !(Y_is_finite)
  )
}

stopifXnotfinite <- function(Y, X, fileConn, msg) {
  result <- this_should_happen(Y, X)
  if (!(result)) {
    # if (!is.null(fileConn)) {
    print_and_append_to_log(c("X is not finite", "\n"), fileConn)
    print_and_append_to_log(c(msg, "\n"), fileConn)
    # }
    stop(msg)
  }
}

# stopifnot(this_should_happen)  # Jul 2019

# custom_sampling_uniform(ir, lb, ub, muj, sdj) {
#   runif(n = length(ir), min = pnorm(q = lb, mean = muj[ir], sd = sdj), max = pnorm(q = ub, mean = muj[ir], sd = sdj))
# }

# quantile_custom_sampling_uniform <- function(ir, lb, ub, muj, sdj) {
#   qnorm(p = runif(n = length(ir), min = pnorm(q = lb, mean = muj[ir], sd = sdj), max = pnorm(q = ub, mean = muj[ir], sd = sdj)), mean = muj[ir], sd = sdj)
# }
#


my_inferCopulaFactorModel <- function(Y, Lambda = diag(ncol(Y)), trueSigma = NULL, nsamp = 1000,
                                      # rand.start = F, 
                                      odens = max(1, round(nsamp / 1000)), impute = any(is.na(Y)),
                                      plugin.threshold = 20,
                                      plugin.marginal = (apply(Y, 2, function(x) {
                                        length(unique(x)) }) > plugin.threshold),
                                      verb = TRUE,
                                      tol = .Machine$double.eps,
                                      first_random_seed = 1000,
                                      random_seed_update_parameter = 10,
                                      output_intermediate = "",
                                      run_id = 0,
                                      fileConn = NULL  # file("output.txt", 'a')
) {
  #
  # This is the main function to perform inference for Gaussian copula factor models.
  # 
  # Args:
  #   Y, a n by p data matrix
  #   Lambda, a p by k matrix (nr variables x nr factors) representing the mapping from factors (columns) to observed variables (rows)
  #   nsamp, No. of samples
  #   odens = for Gibbs-sampling thinning, only Corr matrix after odens sampling is stored
  #   output_intermediate = path used for generating path where intermediate samples covariance are saved
  #   run_id = used for output_intermediate file saving
  #
  #   For details about the arguments, refer to function 'sbgcop.mcmc' in R package 'sbgcop'
  #
  # Return: a list that contains samples of the correlation matrix over all variables (factors + observed variables)
  #   

  print_and_append_to_log(c("Tollerance: ", tol, "\n"), fileConn)  # FG

  # I stores in C.psamp only some of the Gibbs sampled Correlation matrices. 
  # However, in this way to throw-away the burn-in samples from the C.psamp given as output by the function we must throw away the first floor((throw_away_n + 1)/odens_n)
  require(BDgraph)
  library(mvtnorm)
  library(sbgcop)

  random_seed_n <- first_random_seed  # FG for controlling random sampling
  # update_parameter_random_seed_n <- 10  # FG for controlling random sampling

  Y <- as.matrix(Y)
  vnames <- colnames(Y)
  colnames(Y) <- vnames
  # sample size
  n <- dim(Y)[1]
  # No. of observed variables
  p <- dim(Y)[2]

  #### handle Lambda and get prior graph G
  original_lambda_m <- Lambda  # FG
  Lambda <- 1 * (Lambda != 0)
  is_pure_factor_model_b <- (sum(rowSums(Lambda) == 1) == p)
  stopGeneric(is_pure_factor_model_b, fileConn, 'Factor model is not pure')
  # No. of factors
  k <- ncol(Lambda)
  # index of factors with a single indicator. FG used later also as index of the respective variables
  singleton_factor_index <- which(colSums(Lambda) == 1)
  # No. of factors with a single indicator 
  k1 <- length(singleton_factor_index)
  # FG my idea stopifnot(singleton_factor_index == 1:k1)
  exist_singleton_factor_b <- (k1 > 0)
  # index of factors with multiple indicators
  poly_factor_index <- which(colSums(Lambda) > 1)
  rel_sign_orig_lambda <- sign(original_lambda_m)[, poly_factor_index]  # FG used later for adjusting factor sign
  # No. of factors with multiple indicators
  k2 <- length(poly_factor_index)  # TODO actually maybe is assuming that the factors with multiple indicators come after the ones with single indicator. They are in X by construction
  all_lambda_col_nonzero_b <- (k == (k1 + k2))
  stopGeneric(all_lambda_col_nonzero_b, fileConn, 'Loading matrix has some zero columns')

  exist_poly_factor_b <- (k2 > 0)
  if (exist_singleton_factor_b && exist_poly_factor_b) {
    condition_careful <- (max(singleton_factor_index) < min(poly_factor_index))
    stopGeneric(
      condition_careful,
      fileConn,
      'Singleton factors are not all in the first column of loading matrix'
    )
  }
  ## Get the prior graph G
  G1 <- matrix(data = 1, nrow = k, ncol = k) - diag(k)
  if (exist_singleton_factor_b) {
    G2 <- t(Lambda[-singleton_factor_index, ])
  } else {
    G2 <- t(Lambda)
  }
  G3 <- matrix(data = 0, nrow = p - k1, ncol = p - k1)
  G <- rbind(cbind(G1, G2), cbind(t(G2), G3))
  G[lower.tri(G)] <- 0
  ## prior parameters for the G-Wishart distribution
  # Degrees of freedom
  n0 <- p + k2 + 1
  # Scale matrix
  S0 <- diag(p + k2) / n0

  #### Initialize Z, eta, and S
  # FG creates an object R with shape of Y and such that R[i, j] = position of Y[i, j] in sort(unique(Y[, j]))
  # FG so basically R is the ranking of row i in the column j
  R <- NULL
  for (j in 1:p) {
    R <- cbind(R, match(Y[, j], sort(unique(Y[, j]))))
  }
  Rlevels <- apply(R, 2, max, na.rm = TRUE)  # FG takes max of each R column, that means counting unique(Y[, j])
  # TODO what if we use ties.method = "random"?
  # FG Ranks get rank of Y and put average rank for ties, so multiple values of Y could have the same rank
  # Initialize Z
  Ranks <- apply(Y, 2, rank, ties.method = "average", na.last = "keep")
  N <- apply(!is.na(Ranks), 2, sum)
  U <- t(t(Ranks) / (N + 1))
  rm(N)
  Z <- qnorm(U)  # FG matrix of size n, (k1 + k2)
  rm(U)
  # FG Z is later made made my a column-concatenation of [Z_1 , Z_2]
  # FG Z_1 is a n * k1 matrix with the factors with single indicators (corresponding to unstructured variables)
  # FG Z_2 is a n * k2 matrix with the factors with multiple indicators
  # # handle categorical variable
  # Z[, ind.cat] = Y[, ind.cat]
  #
  set.seed(random_seed_n)  # FG set random generator seed
  Zfill <- matrix(data = rnorm(n * p), nrow = n, ncol = p)  # FG randomly initialize Z
  random_seed_n <- update_random_seed(random_seed_n, random_seed_update_parameter)  # FG update random seed
  # Z[is.na(Y)] <- Zfill[is.na(Y)]
  is_na_Y <- is.na(Y)
  Z[is_na_Y] <- Zfill[is_na_Y]
  rm(is_na_Y)
  # pseudo data of factors with a single indicator
  Z1 = eta1 = Z[, singleton_factor_index]  # FG to be changed TODO
  # pseudo data of response variables
  if (exist_singleton_factor_b) {
    Z2 <- Z[, -singleton_factor_index]
  } else {
    Z2 <- Z
  }
  # psuedo data of factors with multiple indicators
  set.seed(random_seed_n)  # FG set random generator seed # eta2 <- matrix(rnorm(n*k2), n)
  eta2 <- matrix(data = rnorm(n * k2), nrow = n, ncol = k2)  # FG
  random_seed_n <- update_random_seed(random_seed_n, random_seed_update_parameter)  # FG update random seed
  eta <- cbind(eta1, eta2)

  # FG  X <- cbind(eta1, eta2, Z2)
  # X = matrix n x (k1 + k2 + (p - k1)) = n x (p + k2)
  # eta1 = matrix n x k1 with transformed data of variables whose factors have single indicators (so 1-2-1 match with variable)
  # eta2 = matrix n x k2 with pseudo-data of factors with multiple indicators
  # Z2  = matrix n x (p - k1) with transformed data of variables whose factors have with multiple indicators
  # If k1 = 0 there is no eta1
  # X = matrix n x (0 + k2 + p) = n x (p + k2)
  X <- cbind(eta, Z2)  # FG so X <- cbind(eta1, eta2, Z2), where if k1 !=0  X <- (Z[, singleton_factor_index], rnorm(n * k2), Z[, -singleton_factor_index])
  X_col_n <- ncol(X)
  col_X_correct <- (X_col_n == (p + k2))
  stopGeneric(col_X_correct, fileConn, c('Initialized Z has', X_col_n, 'instead of', p + k2))
  var_with_poly_factor_in_X_col_index <- (k + 1):X_col_n
  all_var_in_X_col_index <- (1:X_col_n)[-poly_factor_index]
  # FG temporary check
  stopifXnotfinite(Y, X, fileConn, '1st')
  S <- cov(X) + S0
  # S <- cov(X)  # Initialize the covariance matrix
  # if (n < p) {
  #   S <- S + S0  # New wrt to Cui code: add scale matrix to avoid S be singular FG
  # }  # FG
  # LAtest Ruifei version
  # if (rand.start){
  #   S <- solve(rwish(S0*n0, n0 + 1))
  # }else{
  #   S <- cov(X) + S0
  # }
  matrix_not_singular_b <- all(eigen(data.matrix(S))$values != 0)
  stopGeneric(matrix_not_singular_b, fileConn, 'Initialized covariance matrix is singular')

  ####
  Y.pmean <- Y
  # Z.pmean <- Z2   # FG this is not stored anymore
  if (impute) {
    Y.pmean <- matrix(data = 0, nrow = n, ncol = p)
  }
  LPC <- NULL
  stored_samples_n <- floor(nsamp / odens)
  C.psamp <- array(dim = c(p + k2, p + k2, stored_samples_n))  # FG This is used in preparation for the Corr matrix computation, it stores the sampled correlation matrices every odens samples
  # p + k2 = # variables  +  # factors with multiple indicators
  Y.imp <- NULL
  Z.imp <- array(dim = c(n, p, stored_samples_n))
  if (impute) {
    Y.imp <- array(dim = c(n, p, stored_samples_n))
  }
  #dimnames(C.psamp) <- list(colnames(Y), colnames(Y), 1:stored_samples_n)

  last_perc_factor_of_skip_save_n <- 0
  last_perc_factor_of_skip_print_n <- 0
  min_odens_skip_print_and_save_n <- 10
  min_perc_skip_print_n <- 1
  min_perc_skip_save_n <- 10

  # fileConn <- file("output.txt")
  #### start of Gibbs sampling scheme
  # This loop could be put in a function of singleton_factor_index, k, X, S, R, Y, Y.pmean, Y.imp,
  print_and_append_to_log(c("Start Gibbs sampling", "\n"), fileConn)
  for (ns in 1:nsamp) {
    # if (ns <= 2) {
    #   print_and_append_to_log(c("Sampling Iteration: ", ns, "\n"), fileConn)  # FG
    # }
    ## sample Z1 (=eta1)
    # FG these are the factors with single indicator (dummy factors), that are identical to the unstructured variables
    # print_and_append_to_log(c("Start Gibbs sampling", "\n"), fileConn)
    # if (exist_singleton_factor_b) {
    #   print_and_append_to_log(c("Start singleton factor sampling", "\n"), fileConn)
    # }
    for (j_s_factor in singleton_factor_index) {
      # if (!(j %in% ind.cat)) {
      # FG Update Sjc, sdj <- n, muj
      tmp_other_factor_index <- (1:k)[-j_s_factor]  # FG k =  k1 + k2 = nr factors. So these are all the factor numeric id except j_s_factor
      tryCatch(
        {
          Sjc <- S[j_s_factor, tmp_other_factor_index] %*% solve(S[tmp_other_factor_index, tmp_other_factor_index], tol = tol)  # New wrt to Cui code: allow tollerance
        },
        error = function (cond) {
          message(cond)
          stop("Sjc computation failed")
        }
      )
      # compute_sdj <- function(S, j_s_factor, Sjc, tmp_other_factor_index) {
      #   sqrt(S[j_s_factor, j_s_factor] - Sjc %*% S[tmp_other_factor_index, j_s_factor])
      # }
      sdj_n <- sqrt(S[j_s_factor, j_s_factor] - Sjc %*% S[tmp_other_factor_index, j_s_factor])
      muj <- X[, tmp_other_factor_index] %*% t(Sjc)
      # FG simply avoid computing it at every iteration - Could make a binary matrix and 2 binary vectors
      na_in_r_j_bv <- is.na(R[, j_s_factor])  # FG simply avoid computing it at every iteration - actually should be !not_na_in_r_j_bv because R should be as before
      not_na_in_r_j_bv <- !na_in_r_j_bv  # FG simply avoid computing it at every iteration - could actually check is not all TRUE before computing stuff above but it should never happen
      any_na_in_r_j_bv <- any(na_in_r_j_bv)
      # not_na_in_r_j_bv <- !is.na(R[, j]) # FG simply avoid computing it at every iteration - could actually check is not all TRUE before computing stuff above but it should never happen
      # TODO make a function and use it also for X[ir, k2 + j] but it is a bit of a challenge
      # because conditions for max and min must be an input for the function
      if (!plugin.marginal[j_s_factor]) {
        for (r in 1:Rlevels[j_s_factor]) {  # FG loop over 1:length(unique(Y[j_s_factor,]))  (r in 1:length(unique(Y[j_s_factor,])))
          # condition_on_rows_bv <- (R[, j_s_factor]==r & not_na_in_r_j_bv)
          # if (any(condition_on_rows_bv)) {  # FG actually I think this will always be true somehow because of Rlevels
          # ir <- (1:n)[condition_on_rows_bv]
          # FG ir = position in Y of rows i's s.t. Y[i, j_s_factor] has position r in
          #  unique(Y[j_s_factor,]) [ unique(Y[j_s_factor,]) == r ]
          ir <- (1:n)[R[, j_s_factor] == r & not_na_in_r_j_bv]
          # ir <- (1:n)[R[, j_s_factor] == r & !is.na(R[, j_s_factor])]
          # selected_col = j_s_factor, row_value_lower = r - 1, row_value_upper = r + 1,
          # TODO this is not the source of problems: R is derived from Y but then it is applied to X.
          # TODO Bug could be due to repeated values in X?
          lower_bound <- suppressWarnings(max(X[R[, j_s_factor] == r - 1, j_s_factor], na.rm = TRUE))  # FG is the max of X[, j_s_factor] restricted to rows that for Y[, j_s_factor] precede unique(Y[j_s_factor,]) == r
          upper_bound <- suppressWarnings(min(X[R[, j_s_factor] == r + 1, j_s_factor], na.rm = TRUE))
          # TODO problem: when ub = lb = 1 => tmp_norm = +Inf
          # X[ir, j_s_factor] <- quantile_custom_sampling_uniform(ir, lb, ub, muj, sdj)  # FG
          # TODO understand why here
          tmp_min <- pnorm(q = lower_bound, mean = muj[ir], sd = sdj_n)
          tmp_max <- pnorm(q = upper_bound, mean = muj[ir], sd = sdj_n)

          set.seed(random_seed_n)  # FG set random generator seed for runif
          random_sample_unif <- runif(n = length(ir), min = tmp_min, max = tmp_max)
          random_seed_n <- update_random_seed(random_seed_n, random_seed_update_parameter)  # FG update random seed
          tmp_qnorm <- qnorm( p = random_sample_unif, mean = muj[ir], sd = sdj_n )
          X[ir, j_s_factor] <- tmp_qnorm
          if (!this_should_happen(Y, X)) {

            focus_v <- X[R[, j_s_factor] == r - 1, j_s_factor]
            if (any(!is.finite(focus_v))) {
              showData(focus_v, fileConn, "X[R[, j] == r - 1, j]")
            }
            focus_v <- X[R[, j_s_factor] == r + 1, j_s_factor]
            if (any(!is.finite(focus_v))) {
              showData(focus_v, fileConn, "X[R[, j] == r + 1, j]")
            }
            focus_v <- X[, j_s_factor]
            if (any(!is.finite(focus_v))) {
              showData(focus_v, fileConn, "X[, j]")
            }
            focus_v <- R[, j_s_factor]
            if (any(!is.finite(focus_v))) {
              showData(focus_v, fileConn, "R[, j]")
            }

            showData(tmp_qnorm, fileConn, "tmp_qnorm")
            showData(random_sample_unif, fileConn, "random_sample_unif")
            showData(tmp_min, fileConn, "tmp_min")
            showData(tmp_max, fileConn, "tmp_max")
            showData(lower_cound, fileConn, "lower bound")
            showData(upper_cound, fileConn, "upper bound")

            print_and_append_to_log(c("ir", "\n", sep = ""), fileConn)
            print_and_append_to_log(c(capture.output(ir), "\n", sep = ""), fileConn)

            print_and_append_to_log(c("muj[ir]", "\n", sep = ""), fileConn)
            print_and_append_to_log(c(capture.output(muj[ir]), "\n", sep = ""), fileConn)

            print_and_append_to_log(c("sdj", sdj_n, "\n", sep = ""), fileConn)
            # print_and_append_to_log(c(capture.output(sdj), "\n", sep=""), fileConn)
            print_and_append_to_log(c("random_sample_unif", "\n", sep = ""), fileConn)
            print_and_append_to_log(c(capture.output(random_sample_unif), "\n", sep = ""), fileConn)
            stopifXnotfinite(Y, X, fileConn,
              paste(
                "Infinite values while sampling Z1 (singleton factor)", "j=", toString(j), "r=", toString(r), "max_unif", toString(tmp_max), "min_unif", toString(tmp_min)
              )
            )
          }
          # }
        }
      }
      # na_in_r_j_bv <- is.na(R[, j_s_factor])  # FG simply avoid computing it at every iteration - actually should be !not_na_in_r_j_bv because R should be as before
      # ir <- (1:n)[is.na(R[, j_s_factor])]
      # if (any(na_in_r_j_bv)) {  # FG run only if needed
      if (any(na_in_r_j_bv)) {  # FG run only if needed
        ir <- (1:n)[na_in_r_j_bv]
        # ir <- (1:n)[is.na(R[, j_s_factor])]
        set.seed(random_seed_n)  # FG set random generator seed
        X[ir, j_s_factor] <- rnorm(n = length(ir), mean = muj[ir], sd = sdj_n)
        stopifXnotfinite(Y, X, fileConn, "Infinite values while sampling Z1 (singleton factor, missing values)")
        random_seed_n <- update_random_seed(random_seed_n, random_seed_update_parameter)  # FG update random seed
        # }else{
        #   ir <- (1:n)[is.na(R[, j_s_factor])]
        #   X[ir, j_s_factor] = rnorm(length(ir))
        # }
      }
    }
    Z1 = eta1 = X[, singleton_factor_index]  # FG to be changed

    ## sample Z2
    # FG these are the factors whose variables have multiple indicators
    # FG the columns of Z corresponding to Z2 are the last p - k1
    # FG  if (k1 == 0) tmp_var_with_poly_factor_index <- sample(1:p) else tmp_var_with_poly_factor_index <- sample((1:p)[-singleton_factor_index])
    # FG tmp_var_with_poly_factor_index <- index in X of variables whose factors have multiple indicators
    set.seed(random_seed_n)  # FG set random generator seed
    # FG could just sample(k1+1:p-k1) if  singleton_factor_index is (1:k1) like it is indirectly implied
    # FG WARNING instead of using tmp_var_with_poly_factor_index for looping and then always use k2 + j_s_factor, better directly create a set of indices like k2+j_s_factor and loop on it. And then get j_s_factor as Lambda column index somehow
    if (exist_singleton_factor_b) {
      tmp_var_with_poly_factor_index <- sample((1:p)[-singleton_factor_index])  # FG randomly permuted index of variables whose factors have multiple indicators
    } else {
      tmp_var_with_poly_factor_index <- sample(1:p)  # FG randomly permuted 1:p
    }
    random_seed_n <- update_random_seed(random_seed_n, random_seed_update_parameter)  # FG update random seed
    # FG what does k2 + j correspond to?
    # TODO check this
    # FG WARNING j_s_factor is in tmp_var_with_poly_factor_index so is a row index of Lambda (which has size p * k)
    # FG however it is used also as a column index of X and S
    # FG I think this can be done only if singleton_factor_index preceed poly_factor_index
    # FG so that the first k = (k1 + k2) columns of both X and Lambda (and S) are the same
    # if (exist_poly_factor_b) {
    #   print_and_append_to_log(c("Start poly factor sampling", "\n"), fileConn)
    # }
    if (exist_poly_factor_b) {
      for (j_var in tmp_var_with_poly_factor_index) {  # FG Loop over variables that don't have dummy factor assocaited
        # TODO create list which_nonzero_loadings_ls[[j]] <- which(Lambda[j, ] != 0)
        q <- which(Lambda[j_var,] != 0)  # FG j is used as position of Lambda row
        # maybe assert len(q) = 1
        # q = factor corresponding to j
        a <- S[k2 + j_var, q] / S[q, q]  # FG j is used as position of S column
        sdj_n <- sqrt(S[k2 + j_var, k2 + j_var] - a * S[q, k2 + j_var])
        muj <- X[, q] * a
        not_na_in_r_j_bv <- !is.na(R[, j_var]) # FG could actually check is not all TRUE befpre computing stuff above but it should never happen
        if (!plugin.marginal[j_var]) {
          for (r in sort(unique(R[, j_var]))) {
            # ir <- (1:n)[R[, j_var]==r & !is.na(R[, j_var])]
            # condition_on_rows_bv <- (R[, j_var]==r & not_na_in_r_j_bv)
            # if (any(condition_on_rows_bv)) {  # FG actually I think this will always be true somehow because of Rlevels
            # ir <- (1:n)[condition_on_rows_bv]
            ir <- (1:n)[R[, j_var] == r & not_na_in_r_j_bv]
            # selected_col = j, row_value_lower = r - 1, row_value_upper = r + 1,
            lower_bound <- suppressWarnings(max(X[R[, j_var] < r, k2 + j_var], na.rm = TRUE))  # FG j is used as position of X column
            upper_bound <- suppressWarnings(min(X[R[, j_var] > r, k2 + j_var], na.rm = TRUE))  # FG j+k2 , T
            set.seed(random_seed_n)  # FG set random generator seed
            # X[ir, k2 + j] <- quantile_custom_sampling_uniform(ir, lb, ub, muj, sdj)  # FG
            X[ir, k2 + j_var] <- qnorm(
              p = runif(
                n = length(ir),
                min = pnorm(
                  q = lower_bound,
                  mean = muj[ir],
                  sd = sdj_n
                ),
                max = pnorm(
                  q = upper_bound,
                  mean = muj[ir],
                  sd = sdj_n
                )
              ),
              mean = muj[ir],
              sd = sdj_n
            )  # FG j+k2
            stopifXnotfinite(Y, X, fileConn, "Infinite values while sampling Z2 (variables of poly factor)")
            random_seed_n <- update_random_seed(random_seed_n, random_seed_update_parameter)  # FG update random seed
            # }
          }
        }
        # FG simply avoid computing it at every iteration -
        na_in_r_j_bv <- is.na(R[, j_var])  # FG - actually should be !not_na_in_r_j_bv because R should be as before
        # ir <- (1:n)[is.na(R[, j])]
        if (any(na_in_r_j_bv)) {  # FG run only if needed
          ir <- (1:n)[na_in_r_j_bv]
          set.seed(random_seed_n)  # FG set random generator seed
          X[ir, k2 + j_var] <- rnorm(length(ir), muj[ir], sdj_n)  # FG j+k2
          stopifXnotfinite(Y, X, fileConn, "Infinite values while sampling Z2 (variables of poly factor, missing values)")
          random_seed_n <- update_random_seed(random_seed_n, random_seed_update_parameter)  # FG update random seed
        }
      }
      Z2 <- X[, var_with_poly_factor_in_X_col_index]

      ## sample eta2
      A <- S[poly_factor_index, all_var_in_X_col_index] %*% solve(S[all_var_in_X_col_index, all_var_in_X_col_index], tol = tol)  # FG
      sdj_n <- S[poly_factor_index, poly_factor_index] - A %*% S[all_var_in_X_col_index, poly_factor_index]
      # message("sdj size ", dim(sdj))  # FG
      # if (!isSymmetric(sdj)) {
      #   # FG - measure a symmetry of sdj
      #   sdj_transpose <- t(sdj)
      #   sdj_sym <- 0.5 * (sdj + sdj_transpose)
      #   sdj_asym <- 0.5 * (sdj - sdj_transpose)
      #   # alternative: define half_sdj and half_sdj_transpose
      #   sdj_sym_norm <- norm(sdj_sym, "1")
      #   sdj_asym_norm <- norm(sdj_asym, "1")
      #   sym_measure_sdj <- (sdj_sym_norm - sdj_asym_norm) / (sdj_sym_norm + sdj_asym_norm)
      #   # Then −1≤s ≤+1 with the lower bound saturated for an antisymmetric matrox, upper bound saturated for a symmetric one.
      #   message("sdj symmetry measure: ", sym_measure_sdj)
      #   # stopifnot(sym_measure_sdj > 10/100)  # very raw stopping criterion
      #   sdj <- sdj_sym
      # }
      # FG
      muj <- X[, all_var_in_X_col_index] %*% t(A)
      set.seed(random_seed_n)  # FG set random generator seed
      X[, poly_factor_index] <- muj + rmvnorm(n, sigma = sdj_n)  # Generates A Vector Of Random Values From A Multivariate Normal Distribution - not sure because it mentions a dae package
      # stopifXnotfinite(Y, X, fileConn, "6th")
      random_seed_n <- update_random_seed(random_seed_n, random_seed_update_parameter)  # FG update random seed
      # eta2 <- X[, poly_factor_index]  # FG We actually don't need to set this variable
      if (!all(is.finite(X[, poly_factor_index]))) {
        print_and_append_to_log(c("Factors non-singleton:", k2, "\n"), fileConn)
        print_and_append_to_log(c("Input Dataframe", "\n"), fileConn)
        print_and_append_to_log(c(capture.output(Y), "\n", sep = ""), fileConn)
        print_and_append_to_log(c("Current ", "\n"), fileConn)
        print_and_append_to_log(c(capture.output(X[, poly_factor_index]), "\n", sep = ""), fileConn)
      }
      stopifXnotfinite(Y, X, fileConn, "Infinite values while sampling eta2")
    }

    #
    Z <- cbind(Z1, Z2)
    # eta <- cbind(eta1, eta2)  # FG We actually don't need to set this variable

    ### identification condition
    ## Original

    if (exist_poly_factor_b) {
      # print_and_append_to_log(c("Impose identification condition", "\n"), fileConn)
      for (j_factor in poly_factor_index) {
        # FG which(Lambda[, j_factor] != 0)[1] returns a variable, that hence is a number v s.t. 1<=v<=p
        # FG  since the k1 singleton variable/factors are the first k1 columns of X, we must have that v>=k1
        # FG  this implies that k2 + v >= k1 + k2 = k and therefore it corresponds to a column variable of X
        #selected_variable
        concordance_variable <- which(Lambda[, j_factor] != 0)[1]
        is_var <- (concordance_variable > k1)
        stopGeneric(
          is_var,
          fileConn,
          c("Factor", j_factor, "correlation adjustment selected a variable with singleton factor")
        )
        sign_n <- sign(cov(X[, j_factor], X[, k2 + concordance_variable]))
        if (sign_n == 0) {
          msg_v <-c("Factor", j_factor, "correlation sign adjustment set it to zero")
          print_and_append_to_log(c(msg_v, "\n"), fileConn)
          #stopifnot((sign_n != 0))
          stop(cat(msg_v))
        }
        X[, j_factor] <- X[, j_factor] * sign_n  #sign(cov(X[, j], X[, k2 + which(Lambda[, j] != 0)[1]]))
      }
    }
    # Method input: X, k2, poly_factor_index, Lambda, rel_sign_orig_Lambda)
    ## New - written on  20190703
    # # Difference wrt Cui code: now force copula factors to have same signs of input factors
    # # sign( sum(lapply( Lambda[, j] * ,sign)
    # # cov(X[, j], X[, k2 + which(Lambda[, j] != 0)]))
    #
    # # So X[, j] is my \eta_j as long as j in poly_factor_index
    # # Lambda[i, j] = \gamma_{i, j} is the original factor loading for var i and input factor j
    # # Not sure if it should be extended also to singleton_factor_index that is index of both singleton factors and singleton var
    # # sign_cX <- sign(cov(X))  # I am also computing cov btw poly_factor_index stuff that I won't need
    # # rel_sign_cX <- sign_cX[poly_factor_index, -poly_factor_index]
    # #
    # if (k2 > 0) {
    #   rel_sign_cX <- sign(cov(X))[poly_factor_index, -poly_factor_index]  # I just need Cov btw the multi-variable factors and the variables
    #   # # We want sign(diag(sign_cX %*% rel_sign_orig_Lambda))
    #   # # Because we multiply row j of sign_cX with column j of rel_sign_orig_Lambda
    #   # # But this is: for each index
    #   adjust_v <- sign(colSums(t(rel_sign_cX) * rel_sign_orig_lambda))
    #   # # Could this have some zero elements? In theory if the number of variables for a factor is even, it could be
    #   # # Then it is better to edit adjust_v s.t. if it's 0 it is replaced by the original copula sign( cov( X[, j],  X[, k2 + which(Lambda[, j] != 0)[1]] ) )
    #   zero_pos_v <- which(adjust_v == 0)
    #   if (length(zero_pos_v) > 0) {  # Raise something
    #     print_and_append_to_log(
    #       c("There are ", length(zero_pos_v), "factors with sign_adjustment_factor = 0. Force first non-zero factor loading of each multivar factor to be positive.", "\n"),
    #       fileConn
    #     )
    #     # for (jj in zero_pos_v) {
    #     #   adjust_v[jj] <- sign(
    #     #     cov(X[, poly_factor_index[jj]], X[, k2 + which(Lambda[, poly_factor_index[jj]] != 0)[1]])
    #     #   )
    #     # }
    #     adjust_v[zero_pos_v] <- sign(
    #       cov(X[, poly_factor_index[zero_pos_v]], X[, k2 + which(Lambda[, poly_factor_index[zero_pos_v]] != 0)[1]])
    #     )
    #   }
    #   X[, poly_factor_index] <- t(t(X[, poly_factor_index]) * adjust_v)  # multiply each column X[, j] with adjust_v[j]
    # }
    #
    #
    # Could do some elementwise multiplication and then some rows/columns
    # for (jj in 1:length(poly_factor_index){
    #   X[, poly_factor_index[jj]] <- X[, poly_factor_index[jj]] * sign(dot(rel_sign_orig_Lambda[, poly_factor_index[jj]], sign_cX[jj, ]))
    # }

    # 2019
    #
    # FG adjust columns of X corresponding to eta2, looks to me that requires single-identificator factors to come before the others in Lambda columns
    # FG otherwise, change which(Lambda[, j] != 0) in something like
    # for (iPos in (k1+1:k1+k2)) {
    # }
    #
    # TODO double check this j
    # FG WARNING even here j is both an index of Lambda and of X
    # This exploits the fact that the indices j corresponds to those central k2 and hence higher than k1
    # FG  X <- cbind(eta1, eta2, Z2)
    # X = matrix n x (k1 + k2 + (p-k1)) = n x p+k2
    # eta1 = matrix n x k1 with pseudo-data of factors with single indicators (so 1-2-1 with variable)
    # eta2 = matrix n x k2 with pseudo-data of factors with multiple indicators
    # Z2  = matrix n x (p - k1) with transformed data of variables whose factors have with multiple indicators
    # So the Z of Cui_CopulaFactorModel_2018_latest_version.pdf Algorithm 1 is X[, -poly_factor_index] = X[, cat(1:k1, k2+k1+1:p+k2)]
    # if k1=0 => X = matrix n x (0 + k2 + p) = n x p+k2
    #


    # Relocate the mean to zero
    X <- t((t(X) - apply(X, 2, mean)))

    stopifnot(all(is.finite(G)))
    stopifnot(all(is.finite(S0)))
    # stopifnot(all(is.finite(crossprod(X))))
    # crossprod_X <- crossprod(X)
    finite_crossprod_X <- is.finite(crossprod(X))
    stopifnot_condition_X <- all(finite_crossprod_X)

    if (!stopifnot_condition_X) {
      # print_and_append_to_log(c("First random seed:", first_random_seed, "\n"), fileConn)
      print_and_append_to_log(c("Input Dataframe", "\n"), fileConn)
      # print_and_append_to_log(paste(capture.output(Y_df), "\n", sep=""), fileConn)
      print_and_append_to_log(c(capture.output(Y), "\n", sep = ""), fileConn)
      print_and_append_to_log(c("Current X", "\n"), fileConn)
      print_and_append_to_log(c(capture.output(X), "\n", sep = ""), fileConn)
      print_and_append_to_log(c("Current finite_crossprod_X", "\n"), fileConn)
      print_and_append_to_log(c(capture.output(finite_crossprod_X), "\n", sep = ""), fileConn)
      # print_and_append_to_log(toString.finite_crossprod_X), fileConn)
      stop("Some cross-product of X are not finite")
    }
    #stopifnot(stopifnot_condition_X)

    # Sample S - correlation matrix
    set.seed(random_seed_n)  # FG set random generator seed
    if (version_at_least("BDgraph", "2.56")) {
      P <- rgwish(n = 1, adj = G, b = n + n0, D = S0 * n0 + crossprod(X))
    } else if (version_at_least("BDgraph", "2.46")) {
      P <- rgwish(n = 1, adj.g = G, b = n + n0, D = S0 * n0 + crossprod(X))
    } else {
      P <- rgwish(n = 1, adj.g = G, b = n + n0, D = S0 * n0 + crossprod(X))[, , 1]
    }

    random_seed_n <- update_random_seed(random_seed_n, random_seed_update_parameter)  # FG update random seed
    S <- solve(P, tol = tol)  # FG tollerance adde
    S <- cov2cor(S)  # FG transform it to correlation matrix

    # Store the samples every #odens samples
    if (ns %% odens == 0) {
      C <- S #/(sqrt(diag(S)) %*% t(sqrt(diag(S))))
      if (is.null(trueSigma)) {
        lpc <- ldmvnorm(scale(X), C)  # FG compute log of the multivariate normal density
      } else {
        lpc <- ldmvnorm(scale(X), trueSigma)
      }
      LPC <- c(LPC, lpc)
      storing_position_n <- as.integer(ns / odens)

      C.psamp[, , storing_position_n] <- C #cor(X)#
      if (impute) {  # imputing missing values of Y
        # FG missing values of Y are imputed each time
        #
        Y.imp.s <- Y
        for (j in 1:p) {
          # message("Impute ", j, "th column")
          # FG WARNING j is used as an index of both Y and Z columns, but this should be fine
          # FG split up code and execute only when needed
          tmp_isna_Y_bv <- is.na(Y[, j])
          if (any(tmp_isna_Y_bv)) {
            tmp_distro <- pnorm(q = Z[tmp_isna_Y_bv, j], mean = 0, sd = sd(Z[, j]))
            Y.imp.s[tmp_isna_Y_bv, j] <- quantile(
              x = Y[, j],
              probs = tmp_distro,
              na.rm = TRUE,
              names = TRUE,  # FG actually default
              type = 1
            )
          }
          rm(tmp_isna_Y_bv)
          # FG Original version
          # Y.imp.s[is.na(Y[, j]), j] <- quantile(
          #   Y[, j],
          #   pnorm(q = Z[is.na(Y[, j]), j], mean = 0, sd = sd(Z[, j])),
          #   na.rm = TRUE,
          #   type = 1
          # )
        }
        Y.imp[, , storing_position_n] <- Y.imp.s
        # print(Y.imp.s)
        Y.pmean <- ((storing_position_n - 1) / storing_position_n) * Y.pmean + (1 / storing_position_n) * Y.imp.s  # here was the error
        # Y.pmean <- 1 * Y.pmean + storing_position_n ^ (-1) * (Y.imp.s  - Y.pmean)
      }
      Z.imp[, , storing_position_n] <- cbind(Z1, Z2)
    }
    every_x_odens_b <- ((ns %% (odens * min_odens_skip_print_and_save_n)) == 0)
    if (every_x_odens_b) {  # print & and save at most every 10 thinning samples
      perc_completed_n <- round(100 * ns / nsamp, digits = 2)
      perc_factor_of_skip_print_n <- (perc_completed_n %/% min_perc_skip_print_n)  #, digits = 0)
      perc_factor_of_skip_save_n <- (perc_completed_n %/% min_perc_skip_save_n)
    }
    if (
        every_x_odens_b &&
        (verb == TRUE) &&
        (perc_factor_of_skip_print_n > last_perc_factor_of_skip_print_n)
      ) {  # print this every 10 thinning samples
      print_and_append_to_log(
        c(perc_completed_n, "percent done ", date(), "\n"), fileConn
      )
      last_perc_factor_of_skip_print_n <- perc_factor_of_skip_print_n
    }

    if (
        every_x_odens_b &&
        (perc_factor_of_skip_save_n > last_perc_factor_of_skip_save_n) &&
        !stri_isempty(output_intermediate)
      ) {
      saveRDS(
        C.psamp,
        file = sub(
          '.rds',
          paste0('_run_', run_id, '_only_cov_perc_', min_perc_skip_save_n * perc_factor_of_skip_save_n, '_thin_', storing_position_n, '_of_', stored_samples_n, '.rds'),
          output_intermediate
        )
      )
      last_perc_factor_of_skip_save_n <- perc_factor_of_skip_save_n
    }
  }
  # close(fileConn)

  # Z.pmean <- apply(Z.imp, c(1, 2), mean)  # FG this is not stored anymore

  # # 
  #G.ps <- list(Sigma.psamp = C.psamp, Y.pmean = Y.pmean, Y.impute = Y.imp, Z.pmean = Z.pmean, Z.impute = Z.imp, LPC = LPC)
  G.ps <- list(Sigma.psamp = C.psamp)
  class(G.ps) <- "psgc"
  return(G.ps)
}


infer_covariance_and_graph <- function(
  data_df, factor_n, throw_away_odens_n,
  infer_copula_param_ls,
  pc_parameters_ls,
  output_path_pc_algo_obj,
  run2pcalgo_ls,
  run2pcalgo_bad_ls,
  run2suffstat_ls,
  causal_discovery_observation_n = 0,
  iRun_n = 1,
  causal_discovery_algorithm_run_n = 1,
  perform_bootstrap_b = FALSE,
  bootstrap_random_seed_n = 1000,
  bootstrap_random_seed_update_parameter_n = 10,
  fileConn = NULL  # file("output.txt", 'a')
) {
  no_perform_bootstrap_b <- !perform_bootstrap_b
  if (no_perform_bootstrap_b) {
    # cat("No bootstrap performed \n")
    print_and_append_to_log(c("No bootstrap performed \n"), fileConn)
    tmp_run_data_df <- data_df
  } else {
    set.seed(bootstrap_random_seed_n) # for testing set seed to
    print_and_append_to_log(c("Current bootstrap sample:", iRun_n, "of", causal_discovery_algorithm_run_n, "\n"), fileConn)
    print_and_append_to_log(c("Random seed bootstrap sample:", bootstrap_random_seed_n, "\n"), fileConn)
    iteration_rows_id <- sample(nrow(data_df), replace = TRUE)
    bootstrap_random_seed_n <- update_random_seed(bootstrap_random_seed_n, bootstrap_random_seed_update_parameter_n)  # FG update random seed
    tmp_run_data_df <- data_df[iteration_rows_id, ]
    # write.csv(
    #   tmp_run_data_df,
    #   file = '/Users/fabio/rpq_230.csv',
    #   row.names=FALSE,
    #   na=""
    # )
  }

  ## 2.2 Perform inference - first seed of gibbs sampling is always the same, but data changes
  print_and_append_to_log("Perform inference\n", fileConn)
  tmp_infer_copula_param_ls <- infer_copula_param_ls
  tmp_infer_copula_param_ls[['Y']] <- data.matrix(tmp_run_data_df)
  tmp_infer_copula_param_ls[['run_id']] <- iRun_n
  if (no_perform_bootstrap_b) {  # Used for testing, but in theory the repeat would be enough
    tmp_cop_fac_obj <- do.call("my_inferCopulaFactorModel", tmp_infer_copula_param_ls)
  } else {
    repeat {
      tmp_cop_fac_obj <- tryCatch(
        tmp_cop_fac_obj <- do.call("my_inferCopulaFactorModel", tmp_infer_copula_param_ls),
        error = identity
      )
      if (!is(tmp_cop_fac_obj, "error") || no_perform_bootstrap_b) {
        break
      }
      if (is(tmp_cop_fac_obj, "error")) {
        print_and_append_to_log(c("Error:", tmp_cop_fac_obj[[1]], "\n"), fileConn)
        traceback()
        #print_and_append_to_log(c("Error", "\n"), fileConn)
        #for (j_el in tmp_cop_fac_obj) {
        #  if (is.list(j_el)) {
        #    print(j_el)
        #    #out_txt <- paste0(unlist(j_el), collapse = "\n")
        #    #for (k_el in j_el) {
        #    print_and_append_to_log(c(out_txt, "\n"), fileConn)
        #    #}
        #  } else {
        #    print_and_append_to_log(c(j_el, "\n"), fileConn)
        #  }
        #}
      }
      set.seed(bootstrap_random_seed_n) # for testing set seed to
      print_and_append_to_log(c("Repeat current bootstrap sample:", iRun_n, "of", causal_discovery_algorithm_run_n, "\n"), fileConn)
      print_and_append_to_log(c("Random seed bootstrap sample:", bootstrap_random_seed_n, "\n"), fileConn)
      iteration_rows_id <- sample(nrow(data_df), replace = TRUE)
      bootstrap_random_seed_n <- update_random_seed(bootstrap_random_seed_n, bootstrap_random_seed_update_parameter_n)  # FG update random seed
      tmp_run_data_df <- data_df[iteration_rows_id, ]
    }
  }

  # cat(dim(tmp_cop_fac_obj$Sigma.psamp), "\n")
  # Extract samples of the correlation matrix over latent variables, ignoring the first samples (burn-in)
  # C_samples <- tmp_cop_fac_obj$Sigma.psamp[1:factor_n, 1:factor_n, (gibbs_burn_in_n + 1) : gibbs_sampling_n]
  if (throw_away_odens_n < 1) {
    C_samples <- tmp_cop_fac_obj$Sigma.psamp[1:factor_n, 1:factor_n, ]
  } else {
    C_samples <- tmp_cop_fac_obj$Sigma.psamp[1:factor_n, 1:factor_n, -(1:throw_away_odens_n)]
  }
  # the posterior mean
  C <- apply(C_samples, c(1, 2), mean)

  # Allows user to input a fix value to causal_discovery_observation_n
  tmp_causal_discovery_observation_n <- causal_discovery_observation_n
  if (causal_discovery_observation_n == 0) {
    C_sd <- apply(C_samples, c(1, 2), sd)
    # effective sample size
    C_ess <- ((1 - C^2) / C_sd)^2
    # effective_observation_n <- mean(C_ess[upper.tri(C_ess)])
    tmp_causal_discovery_observation_n <- mean(C_ess[upper.tri(C_ess)])
  }

  ## export to CSV
  # write.csv(C, file = args[7])

  # set.seed(random_seed_n)
  #
  #### 2.3 Causal discovery ####
  # MEMO we could directly store the bnlearn version
  ## call the order independent version of the standard PC algorithm

  tmp_suffstat_ls <- list(C = C, n = tmp_causal_discovery_observation_n)
  run2suffstat_ls[[iRun_n]] <- tmp_suffstat_ls
  # For parallelization, we could just make a function returning tmp_suffstat_ls and store it into run2suffstat_ls.
  # After that, in another loop, we can run the causal discovery algorithm and store
  # - run2pcalgo_ls
  # - run2pcalgo_bad_ls

  tmp_pc_parameters_ls <- pc_parameters_ls
  tmp_pc_parameters_ls[['suffStat']] <- tmp_suffstat_ls
  # stopifnot(!identical(old_pc_parameters_ls$suffStat, pc_parameters_ls$suffStat))
  tmp_graph_cfpc <- do.call("pc", tmp_pc_parameters_ls)

  # WARNING we might have to turn off the edge assignment part

  # DEBUG begin
  # Restrict to arc that can be enforced on this graph: they must correspond to egdes that are undirected
  # Check if edge is there
  tmp_blacklist_arcs_absent_b <- TRUE
  if (file.exists(input_path_directed_edges_blacklist)) {
    tmp_amat_cpdag <- as(tmp_graph_cfpc, "amat")
    # DAG/CPDAG (format "cpdag"). Directed egde {from} --> {to} is present <=> amat[{from}, {to}] == 0 AND amat[{to}, {from}] == 1
    tmp_cpdag_blacklist_arcs_present_bv <- (tmp_amat_cpdag[dir_arc_blacklist_mat] == 0) & (tmp_amat_cpdag[dir_arc_blacklist_mat[, rev(colnames(dir_arc_blacklist_mat))]] == 1)
    # MAG/DAG (format "pag"). Directed egde {from} --> {to} is present <=> amat[{from}, {to}] == 2 AND amat[{to}, {from}] == 3
    tmp_pag_blacklist_arcs_present_bv <- (tmp_amat_cpdag[dir_arc_blacklist_mat] == 2) & (tmp_amat_cpdag[dir_arc_blacklist_mat[, rev(colnames(dir_arc_blacklist_mat))]] == 3)
    tmp_blacklist_arcs_present_bv <- (tmp_cpdag_blacklist_arcs_present_bv | tmp_pag_blacklist_arcs_present_bv)
    # tmp_blacklist_arcs_present_bv <- (tmp_amat_cpdag[dir_arc_blacklist_mat] == 1) & (tmp_amat_cpdag[dir_arc_blacklist_mat[, rev(colnames(dir_arc_blacklist_mat))]] == 0)
    tmp_blacklist_arcs_absent_b <- !any(tmp_blacklist_arcs_present_bv)
  }
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
  # bootstrapped_graphnel_list[[iRun_n]] <- as(t(as(with_bgk_amat_cpdag, "matrix")), "graphNEL")  # strange we need to specify matrix

  # In theory we would need to convert to graphNEL only if we add background knowledge, because as.bn would work with both
  # bootstrapped_graphnel_list[[iRun_n]] <- as(t(with_bgk_amat_cpdag), "graphNEL")
  if (tmp_blacklist_arcs_absent_b) {
    run2pcalgo_ls[[iRun_n]] <- tmp_graph_cfpc
  } else {
    run2pcalgo_bad_ls[[iRun_n]] <- tmp_graph_cfpc
    print_and_append_to_log(
      c("Current graph removed because it contained", sum(tmp_blacklist_arcs_present_bv), "blacklisted directed edges", "\n"),
      fileConn
    )
    stopifnot(sum(tmp_blacklist_arcs_present_bv) > 0)
    print_and_append_to_log(dir_arc_blacklist_mat[tmp_blacklist_arcs_present_bv], fileConn)
  }
  # saveRDS(graph_cfpc, file = args[8])
  cat('\n\n')
  # Save temporary output
  if (!stri_isempty(output_path_pc_algo_obj)) {
    tmp_out_path <- paste0(
      file_path_sans_ext(output_path_pc_algo_obj), '_', iRun_n, '.', file_ext(output_path_pc_algo_obj)
    )
    saveRDS(run2pcalgo_ls, file = tmp_out_path)
  }
  return(list(
      bootstrap_random_seed_n = bootstrap_random_seed_n,
      run2suffstat_ls = run2suffstat_ls, run2pcalgo_ls = run2pcalgo_ls, run2pcalgo_bad_ls = run2pcalgo_bad_ls
  ))
}
