# This file defines functions for solving the "mixture distribution
# optimization problem". See the Extreme Deconvolution paper (Bovy,
# Hogg & Roweis, 2011) or the REBayes paper (Koenker & Gu, 2017).

# TO DO: Explain what this function does, and how to use it.
# Refer to Extreme Deconvolution paper for EM algorithm.
mixopt.em <- function (L, w, maxiter = 1e4, tol = 1e-6, verbose = TRUE) {

  # Initialize storage for outputs obj and err.
  obj <- rep(0,maxiter)
  err <- rep(0,maxiter)
    
  # Repeat until convergence criterion is met, or until the maximum
  # number of iterations is reached.
  for (iter in 1:maxiter) {

    # Save the current estimate of the mixture weights.
    w0 <- w

    # E STEP
    
    
    # M STEP

    # CHECK CONVERGENCE
  }

  fit <- list()
  class(fit) <- c("mixopt","list")
  return(fit)
}
