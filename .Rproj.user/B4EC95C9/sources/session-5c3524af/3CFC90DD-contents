#' Estimate the co-occurrence network using the joint approach.
#' @param X Design matrix where the first column is subject ID, second column is index of the time points,
#' third column is measurements of time (if applicable), and remaining columns contain presence/absence
#' of taxa.
#' @param theta.init A matrix containing an initial value for the estimated parameter matrix. If NULL, then
#' this approach will provide an initial estimate for the iterative approach.
#' @param penalty The tuning parameter for the penalized logistic regression. If a number is provided, this
#' value will be used as the penalty. Alternatively, a penalty selection method can be provided. Values allowed
#' are "aic" for AIC, "bic" for BIC, "ebic" for extended BIC, and "mean", "1se", or "max" to use the permutation
#' approach with the mean, one standard error, and maximum, respectively.
#' @param tau.init Initial value for the dampening rate parameter. This must be a numeric.
#' @param tol Value for the tolerance used to determine convergence in the iterative approach. By default, we use 0.0001.
#' @param max.iter Maximum number of iterations in the iterative approach before the algorithm stops. By default, we use 1000.
#' @param nperm If using "mean", "1se", or "max" as the penalty selection method, this is the number of permutations that
#' will be considered. By default, we consider 1000 permutations.
#' @param K If using "mean", "1se", or "max" as the penalty selection method, this is the constant used in front of the
#' selected penalty parameter.
#' @return This function returns a matrix that contains the estimated edges between features in the network. Nonzero
#' elements indicate the presence of an edge, while zero values indicate the absence of an edge.
#' @export
estimate.parameters.j <- function (X, theta.init = NULL, penalty = "1se", tau.init = 0.5, tol = 1e-3, max.iter = 1000,
                                   nperm = 1000, K = 1, symm = "min") {
  N <- nrow(X); p <- ncol(X) - 3; iter <- 1
  X <- as.matrix(X)

  if (is.null(theta.init)) {
    my.vec <- apply(X[,-c(1:3)], MARGIN = 2, FUN = mean)
    theta.init <- logit(my.vec) %>% diag()
    theta.init <- my.vec %>% diag()
  }

  theta.old <- theta.init
  tau.old <- tau.init

  tau.new <- estimate.tau.fixed.theta(X = X, THETA = theta.old)
  theta.new <- estimate.theta.fixed.tau(X = X, tau = tau.new, pen = penalty, nperm = nperm, K = K)

  while (sum(abs(theta.new - theta.old)) > tol && abs(tau.new - tau.old) > tol && iter < max.iter) {
    iter <- iter + 1; theta.old <- theta.new; tau.old <- tau.new

    ### Step 1: Update tau using the value of theta we obtained
    tau.new <- estimate.tau.fixed.theta(X = X, THETA = theta.old)

    ### Step 2: Update theta using the previous value of tau
    theta.new <- estimate.theta.fixed.tau(X = X, tau = tau.new, pen = penalty, nperm = nperm, K = K)
  }

  return(theta.new)
}

