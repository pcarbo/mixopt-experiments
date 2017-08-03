# This file defines functions for solving the "mixture distribution
# optimization problem". See the Extreme Deconvolution paper (Bovy,
# Hogg & Roweis, 2011) or the REBayes paper (Koenker & Gu, 2017).

# Shorthand for machine precision.
eps <- .Machine$double.eps

# Scale each column A[,i] by b[i].
scale.cols <- function (A, b)
  t(t(A) * b)
    
# Compute the mixture distribution objective function given n x k
# conditional likelihood matrix L and mixture weights w, where n is
# the number of samples and k is the number of mixture components.
mixopt.objective <- function (L, w)
  -sum(log(c(L %*% w) + eps))

# TO DO: Explain what this function does, and how to use it.
# Refer to Extreme Deconvolution paper for EM algorithm.
mixopt.em <- function (L, w, maxiter = 1e4, tol = 1e-4, drop.threshold = 1e-8,
                       verbose = TRUE) {

  # Get the number of mixture components.
  k <- ncol(L)
    
  # Initialize the mixture weights.
  if (missing(w))
    w <- rep(1/k,k)
    
  # Initialize storage for outputs obj and err.
  obj <- rep(0,maxiter)
  err <- rep(0,maxiter)

  # Compute the objective function value at the initial iterate.
  f <- mixopt.objective(L,w)
  
  # Repeat until convergence criterion is met, or until the maximum
  # number of iterations is reached.
  if (verbose)
    cat("iter     objective     error\n")
  for (iter in 2:maxiter) {

    # Save the current estimate of the mixture weights and the current
    # objective function value.
    f0 <- f
    w0 <- w

    # E STEP
    # Compute the posterior probabilities
    P <- scale.cols(L,w)
    P <- P / (rowSums(P) + eps)

    # M STEP
    # Update the mixture weights.
    w <- colMeans(P)
    
    # COMPUTE OBJECTIVE
    f <- mixopt.objective(L,w)
    
    # CHECK CONVERGENCE
    # Print the status of the algorithm and check the convergence
    # criterion. Convergence is reached when the maximum difference
    # between the mixture weights at two successive iterations is less
    # than the specified tolerance, or when objective increases.
    err[iter] <- max(abs(w - w0))
    obj[iter] <- f
    if (verbose) {
      progress.str <- sprintf("%4d %+0.6e %0.3e",iter,f,err[iter])
      cat(progress.str)
      cat(rep("\r",nchar(progress.str)))
    }
    if (f > f0) {
      err[iter] <- 0
      obj[iter] <- f0
      w         <- w0
      break
    } else if (err[iter] < tol)
      break
  }
  if (verbose)
    cat("\n")
  
  # Return the fitted model parameters and other optimization info. 
  fit <- list(L = L,w = w,err = err[1:iter],obj = obj[1:iter])
  class(fit) <- c("mixopt","list")
  return(fit)
}

# TO DO: Explain what this function does, and how to use it.
# Refer to REBayes paper for formulation of dual problem.
mixopt.dualip <- function () {

}
