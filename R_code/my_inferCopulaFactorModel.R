update_random_seed <- function(random_seed, parameter = 10) {
  random_seed + parameter
}

version_at_least <- function(pkg, than) {
  as.logical((compareVersion(as.character(packageVersion(pkg)), than) >= 0))
}
# custom_sampling_uniform(ir, lb, ub, muj, sdj) {
#   runif(n = length(ir), min = pnorm(q = lb, mean = muj[ir], sd = sdj), max = pnorm(q = ub, mean = muj[ir], sd = sdj))
# }

# quantile_custom_sampling_uniform(ir, lb, ub, muj, sdj) {
#   qnorm(p = runif(n = length(ir), min = pnorm(q = lb, mean = muj[ir], sd = sdj), max = pnorm(q = ub, mean = muj[ir], sd = sdj)), mean = muj[ir], sd = sdj)
# }
#

my_inferCopulaFactorModel <- function (Y, Lambda = diag(ncol(Y)), trueSigma = NULL, nsamp = 1000, 
                                       odens = max(1, round(nsamp/1000)), impute = any(is.na(Y)), 
                                       plugin.threshold = 20, 
                                       plugin.marginal = (apply(Y, 2, function(x) {
                                       length(unique(x))}) > plugin.threshold), 
                                       verb = TRUE,
                                       tol = .Machine$double.eps,
                                       first_random_seed = 1000,
                                       random_seed_update_parameter = 10,
                                       fileConn = file("output.txt",'a')
                                       ) {
  #
  # This is the main function to perform inference for Gaussian copula factor models.
  # 
  # Args:
  #   Y, a n by p data matrix
  #   Lambda, a p by k matrix representing the mapping from factors to observed variables
  #   nsamp, No. of samples
  #   odens = for Gibbs-sampling thinning, only Corr matrix after odens sampling is stored
  #
  #   For details about the arguments, refer to function 'sbgcop.mcmc' in R package 'sbgcop'
  #
  # Return: a list that contains samples of the correlation matrix over all variables (factors + observed variables)
  #   
  
  
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
  Lambda <- 1 * (Lambda!=0)
  stopifnot(sum(rowSums(Lambda) == 1) == p)  # FG check that is a pure factor model
  # No. of factors
  k <- ncol(Lambda)
  # index of factors with a single indicator
  index.1 <- which(colSums(Lambda) == 1)
  # No. of factors with a single indicator 
  k1 <- length(index.1)
  # FG my idea stopifnot(index.1 == 1:k1)
  # index of factors with multiple indicators
  index.2 <- which(colSums(Lambda) > 1)
  rel_sign_orig_lambda <- sign(original_lambda_m)[, index.2]  # FG used later for adjusting factor sign
  # No. of factors with multiple indicators
  k2 <- length(index.2)  # TODO actually maybe is assuming that the factors with multiple indicators come after the ones with single indicator. They are in X by construction
  ## get the pior graph G
  G1 <- matrix(1, k, k) - diag(k)
  if (k1 == 0) {
    G2 <- t(Lambda)
  } else {
    G2 <- t(Lambda[-index.1, ])
  }
  G3 <- matrix(0, nrow = p-k1, ncol = p-k1)
  G <- rbind(cbind(G1, G2), cbind(t(G2), G3))
  G[lower.tri(G)] <- 0
  ## prior parameters for the G-Wishart distribution
  # degrees of freedom
  n0 <- p + k2 + 1
  # scale matrix
  S0 <- diag(p+k2)/n0
  
  #### initialize Z, eta, and S
  R <- NULL
  for (j in 1:p) {
    R <- cbind(R, match(Y[, j], sort(unique(Y[, j]))))
  }
  Rlevels <- apply(R, 2, max, na.rm = TRUE)
  Ranks <- apply(Y, 2, rank, ties.method = "average", na.last = "keep")
  N <- apply(!is.na(Ranks), 2, sum)
  U <- t(t(Ranks)/(N + 1))
  Z <- qnorm(U)  # FG matrix of size n, (k1 + k2)
  # FG Z is made my a column-concatenation of [Z_1 , Z_2]
  # FG Z_1 is matrix with the factors with single indicators (corresponding to unstructured variables)
  # FG Z_2 is matrix with the factors with multiple indicators
  # # handle categorical variable
  # Z[, ind.cat] = Y[, ind.cat]
  #
  set.seed(random_seed_n)  # FG set random generator seed
  Zfill <- matrix(data = rnorm(n * p), nrow = n, ncol = p)  # FG randomly initialize Z
  random_seed_n <- update_random_seed(random_seed_n, random_seed_update_parameter)  # FG update random seed
  Z[is.na(Y)] <- Zfill[is.na(Y)]
  # pseudo data of factors with a single indicator
  Z1 = eta1 = Z[, index.1]  # FG to be changed TODO
  # pseudo data of response variables
  if (k1 == 0) {
    Z2 <- Z
  } else {
    Z2 <- Z[, -index.1]
  }
  # psuedo data of factors with multiple indicators
  set.seed(random_seed_n)  # FG set random generator seed # eta2 <- matrix(rnorm(n*k2), n)
  eta2 <- matrix(data = rnorm(n * k2), nrow = n, ncol = k2)  # FG
  random_seed_n <- update_random_seed(random_seed_n, random_seed_update_parameter)  # FG update random seed
  eta <- cbind(eta1, eta2)
  X <- cbind(eta, Z2)
  # FG  X <- cbind(eta1, eta2, Z2)
  # X = matrix n x (k1 + k2 + (p-k1)) = n x p+k2
  # eta1 = matrix n x k1 with pseudo-data of factors with single indicators (so 1-2-1 with variable)
  # eta2 = matrix n x k2 with pseudo-data of factors with multiple indicators
  # Z2  = matrix n x (p - k1) with transformed data of variables whose factors have with multiple indicators
  # If k1 = 0 there is no eta1
  # X = matrix n x (0 + k2 + p) = n x p+k2
  # if k1 > 0
  S <- cov(X)  # Initialize the covariance matrix
  if (n < p) {
    S <- S + S0  # New wrt to Cui code: add scale matrix to avoid S be singular FG
  }  # FG
  
  ####
  Y.pmean <- Y
  Z.pmean <- Z2
  if (impute) {
    Y.pmean <- matrix(0, nrow = n, ncol = p)
  }
  LPC <- NULL
  C.psamp <- array(dim = c(p+k2, p+k2, floor(nsamp/odens)))  # FG This is used in preparation for the Corr matrix computation, it stores the sampled correlation matrices evrey odens samples
  # p + k2 = # variables  +  # factors with multiple indicators
  Y.imp <- NULL
  Z.imp <- array(dim = c(n, p, floor(nsamp/odens)))
  if (impute) {
    Y.imp <- array(dim = c(n, p, floor(nsamp/odens)))
  }
  #dimnames(C.psamp) <- list(colnames(Y), colnames(Y), 1:floor(nsamp/odens))

  # fileConn <- file("output.txt")
  #### start of Gibbs sampling scheme
  # This loop could be put in a function of index.1, k, X, S, R, Y, Y.pmean, Y.imp,
  for (ns in 1:nsamp) {
    # message("Sampling Iteration: ", ns)  # FG
    ## sample Z1 (=eta1)
    # FG these are the factors with single indicator (dummy factors), that are identical to the unstructured variables
    for (j in index.1) {
      # if (!(j %in% ind.cat)) {
      ind.tmp <- (1:k)[-j]  # FG k =  k1 + k2 = nr factors
      Sjc <- S[j, ind.tmp] %*% solve(S[ind.tmp, ind.tmp], tol = tol)  # New wrt to Cui code: allow tollerance
      sdj <- sqrt(S[j, j] - Sjc %*% S[ind.tmp, j])
      muj <- X[, ind.tmp] %*% t(Sjc)
      not_na_in_r_j_bv <- !is.na(R[, j]) # FG simply avoid computing it at every iteration - could actually check is not all TRUE before computing stuff above but it should never happen
      # TODO make a function and use it also for X[ir, k2 + j] but it is a bit of a challenge
      # because conditions for max and min must be an input for the function
      if (!plugin.marginal[j]) {
        for (r in 1:Rlevels[j]) {
          # condition_on_rows_bv <- (R[, j]==r & not_na_in_r_j_bv)
          # if (any(condition_on_rows_bv)) {  # FG actually I think this will always be true somehow because of Rlevels
            # ir <- (1:n)[condition_on_rows_bv]
          ir <- (1:n)[R[, j] == r & not_na_in_r_j_bv]
        # ir <- (1:n)[R[, j] == r & !is.na(R[, j])]
          # selected_col = j, row_value_lower = r - 1, row_value_upper = r + 1,
          lb <- suppressWarnings(max(X[R[, j] == r - 1, j], na.rm = TRUE))
          ub <- suppressWarnings(min(X[R[, j] == r + 1, j], na.rm = TRUE))
          set.seed(random_seed_n)  # FG set random generator seed for runif
          # X[ir, j] <- quantile_custom_sampling_uniform(ir, lb, ub, muj, sdj)  # FG
          X[ir, j] <- qnorm(
            p = runif(n = length(ir), min = pnorm(q = lb, mean = muj[ir], sd = sdj), max = pnorm(q = ub, mean = muj[ir], sd = sdj)),
            mean = muj[ir], 
            sd = sdj
          )
          random_seed_n <- update_random_seed(random_seed_n, random_seed_update_parameter)  # FG update random seed
          # }
        }
      }
      na_in_r_j_bv <- is.na(R[, j])  # FG simply avoid computing it at every iteration - actually should be !not_na_in_r_j_bv because R should be as before
      # ir <- (1:n)[is.na(R[, j])]
      if (any(na_in_r_j_bv)) {  # FG run only if needed
        ir <- (1:n)[na_in_r_j_bv]
      # ir <- (1:n)[is.na(R[, j])]
        set.seed(random_seed_n)  # FG set random generator seed
        X[ir, j] <- rnorm(n = length(ir), mean = muj[ir], sd = sdj)
        random_seed_n <- update_random_seed(random_seed_n, random_seed_update_parameter)  # FG update random seed
        # }else{
        #   ir <- (1:n)[is.na(R[, j])]
        #   X[ir, j] = rnorm(length(ir))
        # }
      }
    }
    Z1 = eta1 = X[, index.1]  # FG to be changed
    
    ## sample Z2
    # FG these are the factors whose variables have multiple indicators
    # FG the columns of Z corresponding to Z2 are the last p - k1
    # FG  if (k1 == 0) index.tmp <- sample(1:p) else index.tmp <- sample((1:p)[-index.1])
    # FG index.tmp <- index in X of variables whose factors have multiple indicators
    set.seed(random_seed_n)  # FG set random generator seed
    # FG could just sample(k1+1:p-k1) if  index.1 is (1:k1) like it is indirectly implied
    # FG WARNING instead of using index.tmp for looping and then always use k2 + j, better directly create a set of indices like k2+j and loop on it. And then get j as Lambda column index somehow
    if (k1 == 0) {
      index.tmp <- sample(1:p)  # FG randomly premuted 1:p
    } else {
      index.tmp <- sample((1:p)[-index.1])  # FG randomly permuted index of variables whose factors have multiple indicators
    }
    random_seed_n <- update_random_seed(random_seed_n, random_seed_update_parameter)  # FG update random seed
    # FG what does k2 + j correspond to?
    # TODO check this
    # FG WARNING j is in index.tmp so is a column index of Lambda
    # FG however it is used also as a column index of X and S
    # FG I think this can be done only if index.1 preceed index.2
    # FG so that the first k columns of both X and Lambda (and S) are the same
    if (k2 > 0 ) {
      for (j in index.tmp) {
        q <- which(Lambda[j, ] != 0)  # FG j is used as position of Lambda row
        a <- S[k2+j, q] / S[q, q]  # FG j is used as position of S column
        sdj <- sqrt( S[k2+j, k2+j] - a * S[q, k2+j] )
        muj <- X[, q] * a
        not_na_in_r_j_bv <- !is.na(R[, j]) # FG could actually check is not all TRUE befpre computing stuff above but it should never happen
        if (!plugin.marginal[j]) {
          for (r in sort(unique(R[, j]))) {
            # ir <- (1:n)[R[, j]==r & !is.na(R[, j])]
            # condition_on_rows_bv <- (R[, j]==r & not_na_in_r_j_bv)
            # if (any(condition_on_rows_bv)) {  # FG actually I think this will always be true somehow because of Rlevels
              # ir <- (1:n)[condition_on_rows_bv]
            ir <- (1:n)[R[, j]==r & not_na_in_r_j_bv]
            # selected_col = j, row_value_lower = r - 1, row_value_upper = r + 1,
            lb <- suppressWarnings(max(X[R[, j] < r, k2 + j], na.rm = TRUE))  # FG j is used as position of X column
            ub <- suppressWarnings(min(X[R[, j] > r, k2 + j], na.rm = TRUE))  # FG j+k2 , T
            set.seed(random_seed_n)  # FG set random generator seed
            # X[ir, k2 + j] <- quantile_custom_sampling_uniform(ir, lb, ub, muj, sdj)  # FG
            X[ir, k2 + j] <- qnorm(
              p = runif(n = length(ir), min = pnorm(q = lb, mean = muj[ir], sd = sdj), max = pnorm(q = ub, mean = muj[ir], sd = sdj)),
              mean = muj[ir],
              sd = sdj
            )  # FG j+k2
            random_seed_n <- update_random_seed(random_seed_n, random_seed_update_parameter)  # FG update random seed
            # }
          }
        }
        na_in_r_j_bv <- is.na(R[, j])  # FG - actually should be !not_na_in_r_j_bv because R should be as before
        # ir <- (1:n)[is.na(R[, j])]
        if (any(na_in_r_j_bv)) {  # FG run only if needed
          ir <- (1:n)[na_in_r_j_bv]
          set.seed(random_seed_n)  # FG set random generator seed
          X[ir, k2 + j] <- rnorm(length(ir), muj[ir], sdj)  # FG j+k2
          random_seed_n <- update_random_seed(random_seed_n, random_seed_update_parameter)  # FG update random seed
        }
      }
      Z2 <- X[, (k+1):(k2+p)]
      
      #
      Z <- cbind(Z1, Z2)
      
      ## sample eta2
      ind.tmp <- (1:(p+k2))[-index.2]
      A <- S[index.2, ind.tmp] %*% solve(S[ind.tmp, ind.tmp], tol = tol)  # FG
      sdj <- S[index.2, index.2] - A %*% S[ind.tmp, index.2]
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
      #   # stopifnot(sym_measure_sdj > 10/100)  # very raw stopping driterion
      #   sdj <- sdj_sym
      # }
      # FG
      muj <- X[, ind.tmp] %*% t(A)
      set.seed(random_seed_n)  # FG set random generator seed
      X[, index.2] <- muj + rmvnorm(n, sigma = sdj)  # Generates A Vector Of Random Values From A Multivariate Normal Distribution - not sure because it mentions a dae pacakge
      random_seed_n <- update_random_seed(random_seed_n, random_seed_update_parameter)  # FG update random seed
      eta2 <- X[, index.2]
      
      eta <- cbind(eta1, eta2)

      ### identification condition
      ## Original
      for (j in index.2) {
        X[, j] <- X[, j] * sign(cov(X[, j], X[, k2 + which(Lambda[, j] != 0)[1]]))
      }
      # ## New - written on  20190703
      # # #  difference wrt Cui code: now force copula factors to have same signs of input factors
      # # # sign( sum(lapply( Lambda[, j] * ,sign)
      # # # cov(X[, j], X[, k2 + which(Lambda[, j] != 0)]))
      # #
      # # # So X[, j] is my \eta_j as long as j in index.2
      # # # Lambda[i, j] = \gamma_{i, j} is the original factor loadimg for var i nd input factor j
      # # # Not sure if it should be extended also to index.1 that is index of both singleton factors and signelton var
      # # # sign_cX <- sign(cov(X))  # I am also computing cov btw index.2 stuff that I won't need
      # # # rel_sign_cX <- sign_cX[index.2, -index.2]
      # # #
      # rel_sign_cX <- sign(cov(X))[index.2, -index.2]  # I just need Cov btw the multi-variable factors and the variables
      # # # we want sign(diag(sign_cX %*% rel_sign_orig_Lambda))
      # # # because we multiply row j of sign_cX with column j of rel_sign_orig_Lambda
      # # # But this is: for each index
      # adjust_v <- sign(colSums(t(rel_sign_cX) * rel_sign_orig_lambda))
      # # # Could this have some zero elements? In theory if the number of variables for a factor is even, it could be
      # # # Then it is better to edit adjust_v s.t. if it's 0 it is replaced by the original copula sign( cov( X[, j],  X[, k2 + which(Lambda[, j] != 0)[1]] ) )
      # zero_pos_v = which(adjust_v == 0)
      # if (length(zero_pos_v) > 0) { # TODO raise something
      #   tmp_log_str_v <- c("There are ", length(zero_pos_v), "factors with sign_adjustment_factor = 0. Force first non-zero factor loading of each multivar factor to be positive.", "\n")
      #   cat(tmp_log_str_v)
      #   write(paste(tmp_log_str_v, collapse = " "), file=fileConn, append=TRUE)
      #   for (jj in zero_pos_v) {
      #     adjust_v[jj] <- sign(
      #       cov(X[, index.2[jj]], X[, k2 + which(Lambda[, index.2[jj]] != 0)[1]])
      #     )
      #   }
      # }
      # X[, index.2] <- t(t(X[, index.2]) * adjust_v)  # multiply each column X[, j] with adjust_v[j]
      # # #
      # # # could do some elementwise multiplication and then some rows/columns
      # # # for (jj in 1: length(index.2){
      # # #   X[, index.2[jj]] <- X[, index.2[jj]] * sign( dot(rel_sign_orig_Lambda[, index.2[jj]],  sign_cX[jj, ]))
      # # # }

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
      # So the Z of Cui_CopulaFactorModel_2018_latest_version.pdf Algorithm 1 is X[, -index.2] = X[, cat(1:k1, k2+k1+1:p+k2)]
      # if k1=0 => X = matrix n x (0 + k2 + p) = n x p+k2
      #
    }
    
    # Relocate the mean to zero
    X <- t( (t(X) - apply(X, 2, mean)) )
    
    # Sample S
    set.seed(random_seed_n)  # FG set random generator seed
    stopifnot(all(is.finite(G)))
    stopifnot(all(is.finite(crossprod(X))))
    stopifnot(all(is.finite(S0)))

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
    
    if (ns%%odens == 0) {
      C <- S#/(sqrt(diag(S)) %*% t(sqrt(diag(S))))
      if (is.null(trueSigma)) {
        lpc <- ldmvnorm(scale(X), C)  # FG compute log of the multivariate normal density
      } else {
        lpc <- ldmvnorm(scale(X), trueSigma)
      }
      LPC <- c(LPC, lpc)
      C.psamp[, , ns/odens] <- C #cor(X)#
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
          # FG Original version
          # Y.imp.s[is.na(Y[, j]), j] <- quantile(
          #   Y[, j],
          #   pnorm(q = Z[is.na(Y[, j]), j], mean = 0, sd = sd(Z[, j])),
          #   na.rm = TRUE,
          #   type = 1
          # )
        }
        Y.imp[, , ns/odens] <- Y.imp.s
        # print(Y.imp.s)
        Y.pmean <- ((ns/odens - 1)/(ns/odens)) * Y.pmean + (1/(ns/odens)) * Y.imp.s  # here was the error
      }
      Z.imp[, , ns/odens] <- cbind(Z1, Z2)
    }
    if (verb == TRUE & (ns%%(odens * 10)) == 0) {
      tmp_log_str_v <- c(round(100 * ns/nsamp), "percent done ", date(), "\n")
      cat(tmp_log_str_v)

      write(paste(tmp_log_str_v, collapse = " "), file=fileConn, append=TRUE)
    }
  }
  # close(fileConn)
  
  Z.pmean <- apply(Z.imp, c(1,2), mean)  # Actually we don't use this anymore
  
  # # 
  #G.ps <- list(Sigma.psamp = C.psamp, Y.pmean = Y.pmean, Y.impute = Y.imp, Z.pmean = Z.pmean, Z.impute = Z.imp, LPC = LPC)
  G.ps <- list(Sigma.psamp = C.psamp)
  class(G.ps) <- "psgc"
  return(G.ps)
}
