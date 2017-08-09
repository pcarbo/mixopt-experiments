# This file defines functions for solving the "mixture distribution
# optimization problem". See the Extreme Deconvolution paper (Bovy,
# Hogg & Roweis, 2011) or the REBayes paper (Koenker & Gu, 2017).

# Shorthand for machine precision.
eps <- .Machine$double.eps

# Return the n x n identity matrix.
eye <- function (n)
  diag(rep(1,n))

# Return a sparse n x n identity matrix.
speye <- function (n)
  .symDiagonal(n)

# Return a m x n sparse matrix of all zeros.
spzeros <- function (m, n)
  sparseMatrix(dims = c(m,n),i = NULL,j = NULL)

# Subtract b[i] from each column A[,i].
subtract.cols <- function (A, b)
  t(t(A) - b)

# Scale each column A[,i] by b[i].
scale.cols <- function (A, b)
  t(t(A) * b)
    
# Compute the mixture distribution objective function given n x k
# conditional likelihood matrix L and mixture weights w, where n is
# the number of samples and k is the number of mixture components.
mixopt.objective <- function (L, w)
  -sum(log(drop(L %*% w) + eps))

# TO DO: Explain what this function does, and how to use it.
# Refer to Extreme Deconvolution paper for EM algorithm.
mixopt.em <- function (L, w, maxiter = 1e4, tol = 1e-4, verbose = TRUE) {

  # Get the number of mixture components.
  k <- ncol(L)
    
  # Initialize the mixture weights.
  if (missing(w))
    w <- rep(1/k,k)
    
  # Initialize storage for outputs obj and maxd.
  obj  <- rep(0,maxiter)
  maxd <- rep(0,maxiter)

  # Initialize storage for output "timing".
  timing           <- matrix(0,maxiter,3)
  timing[1,]       <- summary(proc.time())
  colnames(timing) <- names(summary(proc.time()))

  # Compute the objective function value at the initial iterate.
  f <- mixopt.objective(L,w)
  
  # Repeat until convergence criterion is met, or until the maximum
  # number of iterations is reached.
  if (verbose)
    cat("iter     objective max delta\n")
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
    maxd[iter]    <- max(abs(w - w0))
    obj[iter]     <- f
    timing[iter,] <- summary(proc.time())
    if (verbose) {
      progress.str <- sprintf("%4d %+0.6e %0.3e",iter,f,maxd[iter])
      cat(progress.str)
      cat(rep("\r",nchar(progress.str)))
    }
    if (f > f0) {
      maxd[iter] <- 0
      obj[iter]  <- f0
      w          <- w0
      break
    } else if (maxd[iter] < tol)
      break
  }
  if (verbose)
    cat("\n")

  # Reset the timings to zero.
  timing <- timing[1:iter,]
  timing <- subtract.cols(timing,timing[1,])
  
  # Return the fitted model parameters and other optimization info.
  fit <- list(L = L,w = w,maxd = maxd[1:iter],obj = obj[1:iter],
              timing = timing)
  class(fit) <- c("mixopt","list")
  return(fit)
}

# TO DO: Explain what this function does, and how to use it.
# Refer to REBayes paper for formulation of dual problem.
mixopt.dualip <- function (L, maxiter = 1e4, tol = 1e-8, verbose = TRUE) {

  # Get the number of samples.
  n <- nrow(L)
    
  # Get a feasible initial guess for the dual variables.
  x0 <- rep(1/(2*max(L)),n)

  # Solve the dual formulation using the primal-dual interior-point
  # algorithm.
  out <- ipsolver(x = x0,tol = tol,maxiter = maxiter,verbose = verbose,
               
                  # Dual objective.
                  obj = function (x) sum(-log(x + eps)),

                  # Gradient & Hessian of objective.
                  grad = function (x) list(g = -1/(x + eps),
                                           H = spdiag(1/(x^2 + eps))),
                  
                  # Inequality constraints.
                  constr = function (x) c(drop(x %*% L - n),-x),

                  # Jacobian matrix & Hessian of Lagrangian.
                  jac = function (x, z)
                    list(J = rbind(t(L),-speye(n)),
                         W = spzeros(n,n)))
  
  # Recover the dual solution (which gives the mixture weights).
  w <- out$z[1:k]
  
  # Return the fitted model parameters and other optimization info. 
  fit <- list(L = L,w = w,x = out$x,maxd = out$maxd,obj = out$obj,
              ipsolver = out[c("mu","sigma","rx","rc","alpha","ls","x")])
  class(fit) <- c("mixopt","list")
  return(fit)
}
