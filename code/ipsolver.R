# Compute the response of the primal-dual interior-point method merit
# function at (x,z).
ipsolver.merit <- function (x, z, f, d, mu, eps)
  f - sum(d*z) - mu*sum(log(d^2*z + eps))

# Compute the directional derivative of the merit function at (x,z).
ipsolver.gradmerit <- function (x, z, px, pz, g, d, J, mu, eps)
  sum(px * (g - drop(z %*% J - (2*mu/(d - eps)) %*% J))) -
    sum(pz * (d + mu/(z + eps)))

# This function is a simple yet reasonably robust implementation of a
# primal-dual interior-point solver for convex programs with convex
# inequality constraints (it also allows for equality constraints of
# the form Ax = b). It will compute the solution to the following
# optimization problem:
#
#     minimize    f(x)
#     subject to  c(x) < 0,
#
# where f(x) is a convex objective and c(x) is a vector-valued
# function with outputs that are convex in x. There are many
# optimization problems that can be framed in this form (e.g., see the
# book "Convex Optimization" by Boyd & Vandenberghe). The code is
# mostly based on:
#
#    Armand, Gilbert & Jan-Jegou. A feasible BFGS interior point
#    algorithm for solving convex minimization problems. SIAM Journal
#    on Optimization 11, 199-222.
#
# However, to understand what is going on you will need to read up on
# interior-point methods for constrained optimization. A good starting
# point is the "Convex Optmization" book.
#
# Input argument x0 is the initial point for the solver. Input "tol"
# is the tolerance of the convergence criterion; it determines when
# the solver should stop. Input "maxiter" is the maximum number of
# iterations. Input "eqc.tol" is the maximum residual allowable
# residual for the equality constraints Ax = b.
#
# The inputs "obj", "grad", "constr" and "jac" are callback functions
# defined as follows:
# 
#   obj takes input x, the vector of optimization variables, and
#   returns the value of the objective function f(x) at x.
#
#   grad takes the same input x and returns a list with
#   two elements: the gradient and n x n Hessian of the objective at x.
#
#   constr returns the value of the inequality constraint function
#   c(x) at x.
#
#   jac takes two inputs: the primal variables x and the dual
#   variables z. The return value is a list with two elements: the m x
#   n Jacobian matrix (containing the first-order partial derivatives
#   of the inequality constraint functions), and W is the n x n
#   Hessian of the Lagrangian (minus the Hessian of the objective),
#   equal to
#
#     W = z[1]*W1 + z[2]*W2 + ... + z[m]*Wm,
#
#   where Wi is the Hessian of the ith inequality constraint.
#
# There is an optional callback function, "callback", that accepts as
# input the current primal variables (x) and dual variables (z). It is
# called once per iteration and will save all output from this
# function into a list.
#
# If you set "verbose" to true, then at each iteration the solver will
# output the following information: (1) the iteration number; (2)
# objective; (3) barrier parameter, mu; (4) centering parameter,
# sigma; (4) residuals of the perturbed Karush-Kuhn-Tucker system, rx
# & rc; (5) the step size, (6) the number of iterations in the line
# search before a suitable descent step was found, and (7) any output
# provided by the optional callback function.
ipsolver <- function (x, obj, grad, constr, jac, callback = NULL,
                      A = NULL, b = NULL, tol = 1e-8, maxiter = 1e4,
                      newton.solve = "posdef", eqc.tol = 1e-8,
                      verbose = TRUE) {

  # Some algorithm parameters.
  eps      <- 1e-8   # A number close to zero.
  sigmamax <- 0.5    # Maximum centering parameter.
  etamax   <- 0.25   # Maximum forcing number.
  mumin    <- 1e-9   # Minimum barrier parameter.
  alphamax <- 0.995  # Maximum step size.
  alphamin <- 1e-6   # Minimum step size.
  beta     <- 0.75   # Granularity of backtracking search.
  tau      <- 0.01   # Decrease we will accept in line search.
  
  # INITIALIZATION
  # --------------
  # Get the number of primal variables (nv), the number of equality
  # constriants (ne) and the number of inequality constraints (nc).
  nv <- length(x)
  nc <- length(constr(x))
  if (!is.null(A) & !is.null(b)) 
    ne <- length(b)
  else {
    A  <- matrix(0,0,nv)
    b  <- rep(0,0)
    ne <- 0
  }
  
  # Check that x is an interior point and satisfies equality constraints.
  if (any(constr(x) >= 0))
    stop(paste("Initial iterate \"x\" is not an interior point; i.e., does",
               "not satisfy c(x) < 0 for inequality constraints c(x)"))
  else if (ne > 0)
    if (max(abs(A %*% x - b)) > eqc.tol)
      stop(paste("Initial iterate \"x\" does not satisfy equality constraints",
                 "Ax = b within specified tolerance \"eqc.tol\""))
  
  # Initialize the Lagrange multipliers for the inequality constraints
  # (z) and equality constraints (y).
  z <- rep(1,nc)
  y <- rep(0,ne)
  
  # Initialize storage for the outputs.
  out <- list(obj      = rep(0,maxiter),
              maxd     = rep(0,maxiter),
              mu       = rep(0,maxiter),
              sigma    = rep(0,maxiter),
              rx       = rep(0,maxiter),
              rc       = rep(0,maxiter),
              alpha    = rep(0,maxiter),
              ls       = rep(0,maxiter),
              callback = vector("list",maxiter))

  # Initialize storage for the timings output.
  timing           <- matrix(0,maxiter,3)
  timing[1,]       <- summary(proc.time())
  colnames(timing) <- names(summary(proc.time()))
  
  # Repeat while the convergence criterion has not been satisfied, and
  # we haven't reached the maximum number of iterations.
  alpha <- 0
  ls    <- 0
  x0    <- 0
  if (verbose)
    cat(paste("iter       objective  max(d) lg(mu)   sigma  ||rx|| ",
              "||rc|| max(Ax-b) alpha #ls\n"))
  for (iter in 1:maxiter) {

    # Check that x satisfies equality constraints.
    if (ne > 0)
      if (max(abs(A %*% x - b)) > eqc.tol)
        stop(paste("Iterate does not satisfy equality constraints Ax = b",
                   "within specified tolerance \"eqc.tol\""))
  
    # COMPUTE OBJECTIVE, CONSTRAINTS, etc.
    # ------------------------------------
    # Compute the objective (f), the gradient of the objective (g),
    # the Hessian of the objective (H), the inequality constraints (d),
    # the Jacobian of the inequality constraints (J), and the Hessian
    # of the Lagrangian minus the Hessian of the objective (W).
    f <- obj(x)
    d <- constr(x)
    a <- grad(x)
    g <- a$g
    H <- a$H
    a <- jac(x,z)
    J <- a$J
    W <- a$W

    # Compute the unperturbed Karush-Kuhn-Tucker optimality
    # conditions: rx is the dual residual and rc is the
    # complementarity.
    rx <- g + drop(z %*% J + y %*% A)
    rc <- d*z
    r0 <- c(rx,rc)

    # Set some parameters that affect convergence of the primal-dual
    # interior-point method.
    eta        <- min(etamax,norm2(r0)/(nv + nc))
    sigma      <- min(sigmamax,sqrt(norm2(r0)/(nv + nc)))
    dualitygap <- sum(-d*z)
    mu         <- max(mumin,sigma*dualitygap/nc)

    # Print the status of the algorithm.
    maxd <- max(abs(x - x0))
    if (verbose)
      cat(sprintf(paste("%4d %+0.8e %0.1e %+0.3f %0.1e %0.1e %0.1e %0.1e",
                        "%0.1e %03d\n"),
                  iter,f,maxd,log10(mu),sigma,norm2(rx),norm2(rc),
                  max(c(0,abs(A %*% x - b))),alpha,ls))

    # Execute the callback function, if provided.
    if (!is.null(callback))
      out$callback[[iter]] <- callback(x,z)
    
    # Save the status of the algorithm.
    out$obj[iter]   <- f
    out$maxd[iter]  <- maxd
    out$mu[iter]    <- mu
    out$sigma[iter] <- sigma
    out$rx[iter]    <- norm2(rx)
    out$rc[iter]    <- norm2(rc)
    out$alpha[iter] <- alpha
    out$ls[iter]    <- ls
    timing[iter,]   <- summary(proc.time())
    
    # CHECK CONVERGENCE
    # -----------------
    # If the norm of the responses is less than the specified tolerance,
    # we are done. 
    if (norm2(r0)/(nv + nc) < tol)
      break

    # Save the current iterate.
    x0 <- x

    # SOLUTION TO PERTURBED KKT SYSTEM
    # --------------------------------
    if (newton.solve == "posdef") {

      # Compute the search direction of x & z by solving the symmetric
      # positive definite (s.p.d.) linear system Mx = -gb.
      S  <- spdiag(z/(d - eps))
      gb <- g - drop((mu/(d - eps)) %*% J)
      M  <- H + W - t(J) %*% S %*% J
      px <- drop(solve(forceSymmetric(M),-gb))
      pz <- -(z + mu/(b - eps) + drop(S %*% J %*% px))
    } else if (newton.solve == "indef") {

      # Compute the search direction of x & z by solving the
      # indefinite linear system Mx = -r. It is often more efficient
      # (and numerically stable) to solve this system instead of the
      # smaller, s.p.d. system above because the indefinite system is
      # often extremely sparse, whereas the s.p.d. system can be very
      # dense. To view the sparsity pattern of this linear system, use
      # the function from the Matrix package,
      #
      #   image(M)
      #
      # Note that this system can be easily made symmetric by
      # pre-multiplying the linear system by
      #
      #   | I  0 |
      #   | 0 -Z |
      #
      # and likely can be solved much more efficiently using an
      # incomplete Choleksy factorization. However, it seems that this
      # might require a specialized method for solving for saddle
      # point points (see Benzi, Golub & Liesen, 2005, Acta Numerica).
      M <- rbind(cbind(H + W,    t(J)),
                 cbind(-z * J,spdiag(-d + eps)))
      r <- c(rx,-(rc + mu),rep(0,ne))
      if (ne > 0)
        M <- rbind(cbind(M,rbind(t(A),spzeros(nc,ne))),
                   cbind(A,spzeros(ne,nc + ne)))
      p  <- drop(solve(M,-r))
      px <- p[1:nv]
      p  <- p[-(1:nv)]
      pz <- p[1:nc]
      py <- p[-(1:nc)]
    } else
      stop("Input argument \"newton.solve\" should be \"posdef\" or \"indef\"")

    # BACKTRACKING LINE SEARCH
    # ------------------------
    # To ensure global convergence, execute backtracking line search to
    # determine the step length. First, we have to find the largest step
    # size which ensures that z remains feasible. Next, we perform
    # backtracking line search.
    alpha <- alphamax
    i     <- which(z + pz < 0)
    if (length(i) > 0)
      alpha <- alphamax * min(1,min(z[i]/(-pz[i])))

    # Compute the response of the merit function and the directional
    # gradient at the current point and search direction.
    psi  <- ipsolver.merit(x,z,f,d,mu,eps)
    dpsi <- ipsolver.gradmerit(x,z,px,pz,g,d,J,mu,eps)
    ls   <- 0
    while (TRUE) {

      # Compute the candidate iterate, the inequality constraints, and
      # the response of the objective function and merit function at
      # the candidate point.
      ls     <- ls + 1
      xnew   <- x + alpha * px
      ynew   <- y + alpha * py
      znew   <- z + alpha * pz
      f      <- obj(xnew)
      d      <- constr(xnew)
      psinew <- ipsolver.merit(xnew,znew,f,d,mu,eps)

      # Stop backtracking search if we've found a candidate point that
      # sufficiently decreases the merit function and satisfies all the
      # constraints.
      if (sum(d > 0) == 0 & psinew < psi + tau*eta*alpha*dpsi) {
        x <- xnew
        y <- ynew
        z <- znew
        break
      }
      
      # The candidate point does not meet our criteria, so decrease
      # the step size for 0 < beta < 1.
      alpha <- alpha * beta
      if (alpha < alphamin)
        stop("Step size too small")
    }
  }

  # Reset the timings to zero.
  timing <- timing[1:iter,]
  timing <- subtract.cols(timing,timing[1,])
  
  # Return the (primal and dual) solution, and other optimization
  # info.
  return(c(list(x = x,y = y,z = z,timing = timing),
           lapply(out,function (x) x[1:iter])))
}
