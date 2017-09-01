# This file defines functions for solving the "mixture distribution
# optimization problem". See the Extreme Deconvolution paper (Bovy,
# Hogg & Roweis, 2011) or the REBayes paper (Koenker & Gu, 2017).

# Compute the mixture distribution objective function given n x k
# conditional likelihood matrix L and mixture weights w, where n is
# the number of samples and k is the number of mixture components.
mixopt.objective <- function (L, w) {
 if (any(w < 0))
   return(Inf)
 else
   return(-sum(log(drop(L %*% w) + eps)))
}

# Fit a mixture model using EM. Input argument L is the n x k
# conditional likelihood matrix, where n is the number of samples and
# k is the number of mixture components; optional input argument w is
# the initial estimate of the mixture weights.
mixopt.em <- function (L, w, maxiter = 1e4, tol = 1e-4, verbose = TRUE) {

  # Get the number of mixture components.
  k <- ncol(L)
    
  # Initialize the mixture weights.
  if (missing(w))
    w <- rep(1/k,k)
    
  # Initialize storage for outputs obj and maxd.
  obj  <- rep(0,maxiter)
  maxd <- rep(0,maxiter)

  # Initialize storage for the timings output.
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
  class(fit) <- c("mixopt.em","list")
  return(fit)
}

# Fit a mixture model by solving the primal problem using a
# primal-dual interior-point method; see ipsolver.R for algorithm
# details. Input argument L is the n x k conditional likelihood
# matrix, where n is the number of samples and k is the number of
# mixture components
mixopt.ip <- function (L, maxiter = 1e4, tol = 1e-8, verbose = TRUE) {

  # Get the number of mixture components.
  k <- ncol(L)
  n <- nrow(L)
  
  # Solve the dual formulation using the primal-dual interior-point
  # algorithm. Note that the indefinite system for solving the Newton
  # step is very sparse, so we set newton.solve = "indef".
  out <- ipsolver(x = rep(1/k,k),tol = tol,maxiter = maxiter,
                  verbose = verbose,newton.solve = "indef",
                  
                  # Objective.
                  obj = function (x) mixopt.objective(L,x) + n*sum(x),

                  # Gradient & Hessian of objective.
                  grad = function (x) {
                    y <- c(L %*% x)
                    return(list(g = n - colSums(L/(y + eps)),
                                H = diagsq(L,1/y)))
                  },
                  
                  # Inequality constraints.
                  constr = function (x) -x,

                  # Jacobian matrix & Hessian of Lagrangian.
                  jac = function (x, z)
                    list(J = -speye(k),
                         W = spzeros(k,k)))

  # Get the normalized mixture weights.
  w <- out$x
  w <- w/sum(w)
  
  # Return the fitted model parameters and other info returned by the
  # optimization algorithm.
  fit <-
    list(L = L,w = w,obj = out$obj,maxd = out$maxd,timing = out$timing,
         ipsolver = out[setdiff(names(out),c("max","timing","obj"))])
  class(fit) <- c("mixopt.ip","list")
  return(fit)
}

# Fit a mixture model by solving the dual problem using a primal-dual
# interior-point method; see ipsolver.R for algorithm details, and see
# the REBayes paper for formulation of dual problem. Input argument L
# is the n x k conditional likelihood matrix, where n is the number of
# samples and k is the number of mixture components
mixopt.dualip <- function (L, maxiter = 1e4, tol = 1e-8, verbose = TRUE) {

  # Get the number of samples.
  n <- nrow(L)
    
  # Get a feasible initial guess for the dual variables.
  x0 <- rep(1/(2*max(L)),n)

  # Solve the dual formulation using the primal-dual interior-point
  # algorithm. Note that the indefinite system for solving the Newton
  # step is very sparse, so we set newton.solve = "indef".
  out <- ipsolver(x = x0,tol = tol,maxiter = maxiter,
                  verbose = verbose,newton.solve = "indef",
               
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

  # Recover the dual solution (which gives the mixture weights), and
  # return the fitted model parameters and other optimization info
  # returned by the optimization algorithm.
  w   <- out$z[1:k]
  fit <- list(L = L,w = w,obj = mixopt.objective(L,w),timing = out$timing,
              ipsolver = out[setdiff(names(out),"timing")])
  class(fit) <- c("mixopt.dualip","list")
  return(fit)
}

