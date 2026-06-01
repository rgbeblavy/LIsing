library(tidyverse); library(glmnet); library(MASS); library(gtools); library(SIS)

#' Compute the dampening rate parameter for any subject given the times of their repeated measurements.
#' @param time.vec A vector of the observation times (or indices of repeated measurements.)
#' @param tau The value of the dampening rate parameter to compute the dampening rate parameter matrix for.
#' This should be a numeric.
#' @param C The constant to include in the model. By default, we use C = 5.
#' @return This function returns a matrix that contains the dampening rates for a given subject using
#' the exponential model.
#' @export
get.phi <- function (time.vec, tau, C = 5) {
  phi <- matrix(data = 0, nrow = length(time.vec), ncol = length(time.vec))

  for (j in 1:length(time.vec)) {
    for (j.prime in 1:length(time.vec)) {
      phi[j, j.prime] <- exp(-C * tau * (time.vec[j] - time.vec[j.prime])^2)
    }
  }

  return(phi)
}

#' This function finds the position in the parameter vector when transforming the parameter matrix to a vector for all off-diagonal elements.
#' @param s A numeric given the row index of interest in the parameter matrix.
#' @param t A numeric given the column index of interest in the parameter matrix.
#' @return This function returns the index in the parameter vector for off-diagonal elements of the parameter matrix.
#' @export
get.index <- function (s, t) {
  ### s and t are the nodes of the parameter we are interested in
  ### Note: This only works for off-diagonal elements
  j <- min(s,t) + ( ((max(s,t) - 2)*(max(s,t) - 1)) / 2)
  return(j)
}

#' This function converts the parameter vector estimated by our joint approach to a parameter matrix.
#' @param theta.tilda Estimated parameter vector.
#' @param p Number of nodes. This must be a numeric.
#' @return This function returns the index in the parameter matrix given the transformed parameter vector.
#' @export
vec.to.matrix <- function (theta.tilda, p) {
  theta.mat <- matrix(data = 0, nrow = p, ncol = p); m <- choose(p, 2)

  # Step 1: Fill in off-diagonal elements
  for (s in 1:p) {
    for (t in 1:p) {

      if (s != t) {
        j <- get.index(s, t)

        theta.mat[s,t] <- theta.tilda[j]
      }
    }
  }

  # Step 2: Fill in diagonal elements
  diag(theta.mat) <- theta.tilda[(m+1):(m+p)]

  return(theta.mat)
}

#' This function generates the kth standard basis vector in p-dimensional real space.
#' @param p A numeric giving the total length of the vector.
#' @param k A numeric giving the position the 1 should be placed to generate a basis vector.
#' @return This function returns a vector with 1 in position k and 0s elsewhere.
#' @export
make.basis <- function (p, k) {
  replace(numeric(p), k, 1)
}

#' This function transforms the parameter matrix theta into a vector.
#' @param THETA Estimated parameter matrix.
#' @return This function returns a vector containing the same values in the parameter matrix provided.
#' @export
transformed.theta <- function(THETA) {
  p <- dim(THETA)[1]
  m <- choose(p, 2)

  theta.tilda <- vector(mode = "logical", length = m + p)

  ## Fill in off-diagonal elements
  for (s in 1:p) {
    for (t in s:p) {

      if (s != t) {
        j <- get.index(s,t)
        theta.tilda[j] <- THETA[s,t]
      }

    }
  }

  ## Fill in diagonal elements
  for (j in (m+1):(m+p)) {
    theta.tilda[j] <- THETA[j-m,j-m]
  }

  return(theta.tilda)
}

#' This function generates the transformed design matrix for the joint approach given a fixed value
#' for the dampening rate parameter.
#' @param X Original design matrix containing presence/absence data for all subjects across all time points.
#' The first column of X should be "ID" for the subject IDs, the second column of X should be
#' "Measure" to index repeated measurements on each subject, and the third column of X should be
#' "Time" to measure time points for each subject. The remaining columns must have presence/absence data for the features.
#' @param tau Fixed value for the dampening rate parameter. This must be a numeric.
#' @return This function returns the design matrix to be used in the logistic regression problem for our joint approach.
#' @export
transformed.model.matrix <- function (X, tau) {
  ## X is the original design matrix for longitudinal settings
  ## tau is the dampening rate parameter
  ## Column 1 of X should be "ID" for the subject IDs
  ## Column 2 of X should be "Measure" to index repeated measurements on each subject
  ## Column 3 of X should be "Time" to measure time points for each subject
  N <- nrow(X); p <- ncol(X) - 3; n <- length(unique(X[,1]))
  m <- choose(p, 2)

  n.vec <- vector(mode = "logical", length = n)
  Z.tilda <- matrix(data = 0, nrow = N*p, ncol = m + p)

  for (s in 1:p) {

    for (i in 1:n) {

      Xi <- X[which(X[,1] == i),]

      ### The chunk below is for the case where we have just one repeated measurement on a given subject
      if (is.vector(Xi)) {
        n.vec[i] <- 1; ni <- n.vec[i] # if Xi becomes a vector then we only had 1 repeated measure on that subject

        if (i == 1) {
          u <- N*(s-1) + 1
        } else {
          u <- N*(s-1) + sum(n.vec[1:(i-1)]) + 1
        }

        for (t in 1:p) {

          if (s != t) {
            v <- get.index(s,t)

            Z.tilda[u,v] <- Xi[t+3]

          }

        } # end loop over nodes "t"

        for (v in (m+1):(m+p)) {
          Z.tilda[u,v] <- Xi[s+3] * (s == v - m)
        }

        ### The chunk below is for the case where we have multiple repeated measurements on a given subject
      } else {
        n.vec[i] <- Xi[,2] %>% unique() %>% length(); ni <- n.vec[i]

        PHI <- get.phi(Xi[,3], tau = tau)

        for (j in 1:ni) {

          if (i == 1) {
            u <- N*(s-1) + j
          } else {
            u <- N*(s-1) + sum(n.vec[1:(i-1)]) + j
          }

          for (t in 1:p) {

            if (s != t) {
              v <- get.index(s,t)

              Z.tilda[u,v] <- as.numeric(t(make.basis(ni, j)) %*% PHI %*% Xi[,t+3])

            }

          } # end loop over nodes "t"

          for (v in (m+1):(m+p)) {
            my.sum <- 0

            for (j.prime in 1:ni) {
              if (j != j.prime) {
                my.sum <- my.sum + PHI[j,j.prime]*Xi[j.prime,s+3]
              }
            }

            Z.tilda[u,v] <- (1 + my.sum) * (s == v - m)
          }

        } # end loop over repeated measures
      } # end conditional statement


    } # end loop over subjects
  } # end loop over nodes

  return(Z.tilda)
}

#' This function estimates the edge parameter matrix using nodewise regression for a fixed value of the dampening
#' rate parameter.
#' @param X Original design matrix containing presence/absence data for all subjects across all time points.
#' The first column of X should be "ID" for the subject IDs, the second column of X should be
#' "Measure" to index repeated measurements on each subject, and the third column of X should be
#' "Time" to measure time points for each subject. The remaining columns must have presence/absence data for the features.
#' @param tau Fixed value for the dampening rate parameter. This must be a numeric.
#' @param pen The tuning parameter for the penalized logistic regression. If a number is provided, this
#' value will be used as the penalty. Alternatively, a penalty selection method can be provided. Values allowed
#' are "aic" for AIC, "bic" for BIC, "ebic" for extended BIC, and "mean", "1se", or "max" to use the permutation
#' approach with the mean, one standard error, and maximum, respectively.
#' @param nperm If using "mean", "1se", or "max" as the penalty selection method, this is the number of permutations that
#' will be considered. By default, we consider 1000 permutations.
#' @param K If using "mean", "1se", or "max" as the penalty selection method, this is the constant used in front of the
#' selected penalty parameter.
#' @param symm This is the method used to symmetrize the parameter matrix. This can be either "min" or "max" depending
#' on if the maximum or minimum rule is desired.
#' @return This function returns the estimated edge parameters between features for a fixed value of the dampening rate parameter.
#' @export
estimate.theta.nw <- function (X, tau, pen, nperm = 100, K = 1, symm = "min") {
  ### Defining response and model matrix for joint methods
  y.tilda <- as.vector(as.matrix(X[,-c(1:3)]))
  Z.tilda <- transformed.model.matrix(X, tau = tau)

  ### Implement longitudinal nodewise regression
  N <- nrow(X)
  n <- X[,1] %>% unique() %>% length(); p <- ncol(X) - 3; m <- choose(p, 2)
  theta.est <- matrix(data = 0, nrow = p, ncol = p)

  for (s in 1:p) {
    # get transformed response
    y_s <- as.vector(X[,s+3]) # response

    # get row indices of transformed model matrix
    row.idx <- NULL; sum_ni <- 0

    for (i in 1:n) {
      Xi <- X[which(X[,1] == i),]

      if (is.vector(Xi)) {
        n_i <- 1 # only one repeated measurement if above variable is a vector
      } else {
        n_i <- nrow(Xi)
      }

      for (j in 1:n_i) {
        row.idx <- c(row.idx, N*(s-1) + sum_ni + j)
      }

      sum_ni <- sum_ni + n_i

    }

    # get column indices of transformed model matrix
    col.idx <- NULL

    for (t in setdiff(1:p, s)) {
      col.idx <- c(col.idx, get.index(s,t))
    }

    Z_s <- Z.tilda[row.idx,c(col.idx, m+s)]

    # Selecting penalty
    if (is.numeric(pen)) {
      opt.pen <- pen

    } else if (pen %in% c("aic", "bic", "ebic")) {
      opt.pen <- tune.fit(x = Z_s, y = y_s, family = "binomial", penalty = "lasso", tune = pen)$lambda

    } else if (pen %in% c("mean", "1se", "max")) {

      pens <- vector(mode = "logical", length = nperm)

      for (j in 1:nperm) {
        yperm <- sample(y_s, replace = FALSE, size = length(y_s))

        pens[j] <- max(abs( ((t(Z_s) %*% (yperm - (1/2))) / length(yperm)) ))
      }

      if (pen == "mean") {
        opt.pen <- K*mean(pens)
      } else if (pen == "1se") {
        opt.pen <- K*(mean(pens) + sd(pens))
      } else if (pen == "max") {
        opt.pen <- K*max(pens)
      }

    } else {
      return("Error: Invalid penalty criterion provided as input.")
    }

    # Network estimation
    res <- glmnet(x = Z_s, y = y_s, family = "binomial", alpha = 1, lambda = opt.pen, intercept = FALSE)
    est.coefs <- coef(res)[-1] %>% as.vector()

    # Store these estimates in a matrix
    theta.est[s,s] <- est.coefs[length(est.coefs)]
    theta.est[s,-c(s)] <- est.coefs[1:(length(est.coefs)-1)]
  }



  ### Symmetrize parameter matrix
  if (symm == "max") {
    theta.est.wmax <- matrix(data = 0, nrow = p, ncol = p)
    diag(theta.est.wmax) <- diag(theta.est)

    for (s in 1:p) {
      for (t in s:p) {

        if (abs(theta.est[s,t]) > abs(theta.est[t,s])) {
          theta.est.wmax[s,t] <- theta.est[s,t]
          theta.est.wmax[t,s] <- theta.est[s,t]
        } else {
          theta.est.wmax[s,t] <- theta.est[t,s]
          theta.est.wmax[t,s] <- theta.est[t,s]
        }

      }
    }

    return(theta.est.wmax)

  } else if (symm == "min") {
    theta.est.wmin <- matrix(data = 0, nrow = p, ncol = p)
    diag(theta.est.wmin) <- diag(theta.est)

    for (s in 1:p) {
      for (t in s:p) {

        if ( abs(theta.est[s,t]) < abs(theta.est[t,s]) ) {
          theta.est.wmin[s,t] <- theta.est[s,t]
          theta.est.wmin[t,s] <- theta.est[s,t]
        } else {
          theta.est.wmin[s,t] <- theta.est[t,s]
          theta.est.wmin[t,s] <- theta.est[t,s]
        }

      }
    }

    return(theta.est.wmin)

  } else {
    return("Error: Must use max or min rule.")
  }

}

#' This function estimates the edge parameter matrix using the joint approach for a fixed value of the dampening
#' rate parameter.
#' @param X Original design matrix containing presence/absence data for all subjects across all time points.
#' The first column of X should be "ID" for the subject IDs, the second column of X should be
#' "Measure" to index repeated measurements on each subject, and the third column of X should be
#' "Time" to measure time points for each subject. The remaining columns must have presence/absence data for the features.
#' @param tau Fixed value for the dampening rate parameter. This must be a numeric.
#' @param pen The tuning parameter for the penalized logistic regression. If a number is provided, this
#' value will be used as the penalty. Alternatively, a penalty selection method can be provided. Values allowed
#' are "aic" for AIC, "bic" for BIC, "ebic" for extended BIC, and "mean", "1se", or "max" to use the permutation
#' approach with the mean, one standard error, and maximum, respectively.
#' @param nperm If using "mean", "1se", or "max" as the penalty selection method, this is the number of permutations that
#' will be considered. By default, we consider 1000 permutations.
#' @param K If using "mean", "1se", or "max" as the penalty selection method, this is the constant used in front of the
#' selected penalty parameter.
#' @return This function returns the estimated edge parameters between features for a fixed value of the dampening rate parameter.
#' @export
estimate.theta.fixed.tau <- function (X, tau, pen, nperm = 1000, K = 1) {
  N <- nrow(X); n <- X[,1] %>% unique() %>% length(); p <- ncol(X) - 3; m <- choose(p, 2)

  ### Defining response and model matrix for joint methods
  y.tilda <- as.vector(as.matrix(X[,-c(1:3)]))
  Z.tilda <- transformed.model.matrix(X, tau = tau)

  ### Determine optimal penalty
  # Selecting penalty
  if (is.numeric(pen)) {
    opt.pen <- pen

  } else if (pen %in% c("aic", "bic", "ebic")) {
    opt.pen <- tune.fit(x = Z.tilda, y = y.tilda, family = "binomial", penalty = "lasso", tune = pen)$lambda

  } else if (pen %in% c("mean", "1se", "max")) {

    pens <- vector(mode = "logical", length = nperm)

    for (j in 1:nperm) {
      yperm <- sample(y.tilda, replace = FALSE, size = length(y.tilda))

      pens[j] <- max(abs( ((t(Z.tilda) %*% (yperm - (1/2))) / length(yperm)) ))
    }

    if (pen == "mean") {
      opt.pen <- K*mean(pens)
    } else if (pen == "1se") {
      opt.pen <- K*(mean(pens) + sd(pens))
    } else if (pen == "max") {
      opt.pen <- K*max(pens)
    }

  } else {
    return("Error: Invalid penalty criterion provided as input.")
  }

  ### Implement longitudinal nodewise regression
  res <- glmnet(x = Z.tilda, y = y.tilda, family = "binomial", alpha = 1,
                lambda = opt.pen, intercept = FALSE)
  est.coefs <- coef(res)[-1] %>% as.vector()
  theta.est <- vec.to.matrix(est.coefs, p = p)

  return(theta.est)
}

#' This function computes the value of the joint log-likelihood function given the edge parameter matrix, dampening rate
#' parameter, and fixed design matrix.
#' @param X Original design matrix containing presence/absence data for all subjects across all time points.
#' The first column of X should be "ID" for the subject IDs, the second column of X should be
#' "Measure" to index repeated measurements on each subject, and the third column of X should be
#' "Time" to measure time points for each subject. The remaining columns must have presence/absence data for the features.
#' @param THETA Value for the edge parameter matrix. This must be a symmetric matrix where the number of rows and
#' columns are equal to the number of features.
#' @param tau Fixed value for the dampening rate parameter. This must be a numeric.
#' @return This function returns a numeric with the value for the log-likelihood.
#' @export
longitudinal.lik <- function (THETA, tau, X) {
  n <- X[,1] %>% unique() %>% length()
  N <- nrow(X); p <- ncol(X) - 3

  my.sum1 <- 0 # Initialize a counter for the first summation
  my.sum2 <- 0 # Initialize a counter for the second summation

  for (i in 1:n) {
    Xi <- X[which(X[,1] == i),]

    if (is.vector(Xi)) {
      ni <- 1; xi <- as.vector(t(Xi[-c(1:3)])); PHI <- 1

      my.sum1 <- my.sum1 + as.numeric(t(xi) %*% kronecker(PHI, THETA) %*% xi)

      for (s in 1:p) {
        my.num2 <- 0; Xi.tmp <- Xi; Xi.tmp[s+3] <- 1; ind <- which(THETA[s,] != 0)

        for (t in ind) {
          my.num2 <- my.num2 + as.numeric(Xi.tmp[t+3] * THETA[s,t])
        }

        my.sum2 <- my.sum2 + log(1 + exp(my.num2))
      }

    } else {
      ni <- nrow(Xi); xi <- as.vector(t(Xi[,-c(1:3)])); PHI <- get.phi(Xi[,3], tau = tau)

      my.sum1 <- my.sum1 + as.numeric(t(xi) %*% kronecker(PHI, THETA) %*% xi)

      for (j in 1:ni) {
        ej <- make.basis(ni, j)

        for (s in 1:p) {
          my.num2 <- 0; Xi.tmp <- Xi; Xi.tmp[j,s+3] <- 1; ind <- which(THETA[s,] != 0)

          for (t in ind) {
            my.num2 <- my.num2 + as.numeric(t(ej) %*% PHI %*% Xi.tmp[,t+3]) * THETA[s,t]
          }

          my.sum2 <- my.sum2 + log(1 + exp(my.num2))
        }

      }

    }

  }

  my.sum <- my.sum1 - my.sum2

  return(my.sum)
}

#' This function computes the optimal value of the dampening rate parameter given a fixed edge parameter matrix.
#' @param X Original design matrix containing presence/absence data for all subjects across all time points.
#' The first column of X should be "ID" for the subject IDs, the second column of X should be
#' "Measure" to index repeated measurements on each subject, and the third column of X should be
#' "Time" to measure time points for each subject. The remaining columns must have presence/absence data for the features.
#' @param THETA Fixed value for the edge parameter matrix. This must be a symmetric matrix where the number of rows and
#' columns are equal to the number of features.
#' @param s The step size to be used for the bisection approach.
#' #' @param tol Value for the tolerance used to determine convergence in the iterative approach. By default, we use 0.0001.
#' @param max.iter Maximum number of iterations in the iterative approach before the algorithm stops. By default, we use 100.
#' @return This function returns a numeric with the value for the log-likelihood.
#' @export
estimate.tau.fixed.theta <- function(X, THETA, s = 0.1, tol = 1e-3, max.iter = 100) {
  tau.old <- 1; lik.vals <- NULL; iter <- 1

  tau.grid <- c(0.01, seq(0.1, 1, s))

  for (tau in tau.grid) {
    lik.vals <- c(lik.vals,
                  longitudinal.lik(THETA, tau, X))
  }

  tau.new <- tau.grid[which(lik.vals == max(lik.vals))] %>% max()

  ### Iteratively apply bisection approach until convergence
  while (abs(tau.new - tau.old) > tol && iter < max.iter) {
    iter <- iter + 1; tau.old <- tau.new; lik.vals <- NULL

    tau.grid <- seq(tau.old - s, tau.old + s, s/2)
    tau.grid <- tau.grid[which(tau.grid > 0)] # be sure to exclude 0 values

    for (tau in tau.grid) {
      lik.vals <- c(lik.vals,
                    longitudinal.lik(THETA, tau, X))
    }

    tau.new <- tau.grid[which(lik.vals == max(lik.vals))] %>% max()

    s <- s/2
  }

  return(tau.new)
}



