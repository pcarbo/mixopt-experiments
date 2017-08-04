# This short script demonstrates the use of the interior-point solver to
# compute the solution to a quadratic program with convex objective
# (i.e. positive-definite Hessian) and convex, quadratic inequality
# constraints. More precisely, it finds the solution to the following
# optimization problem:
#
#   minimize    x'*H*x/2 + u'*x
#   subject to  ci(x) < b
#
# where the inequality constraints are quadratic functions:
#
#   ci(x) = x'*Pi*x/2 + ri'*x
#
# This particular example originally comes from the book: H. P. Schwefel
# (1995) Evolution and Optimum Seeking. The minimium occurs at (0,1,2,-1).
source("../code/ipsolver.R")

# Return the quadratic (L2) norm of x with respect to matrix A.
qnorm <- function (x, A)
  c(sqrt(sum(x * (A %*% x))))

# Return x'*A*x/2 + b'*x + c.
quadf <- function (x, A, b, c = 0)
  qnorm(x,A)^2/2 + sum(b*x) + c

# Compute the m x n Jacobian of the vector-valued constraint function,
# and the the n x n Hessian of the Lagrangian (minus the Hessian of the
# objective) for the quadratic program described in the comments above.
qpjacobian <- function (x, z, constraints) {
  n <- length(x)
  m <- length(constraints)
  J <- matrix(0,m,n)
  W <- matrix(0,n,n)
  for (i in 1:m) {
    a     <- constraints[[i]]
    J[i,] <- c(a$P %*% x + a$r)
    W     <- W + z[i]*a$P
  }
  return(list(J = J,W = W))
}

# Define the quadratic objective function.
H <- diag(c(2,2,4,2))
u <- c(-5,-5,-21,7)

# Define the quadratic inequality constraints.
constraints <-
  list(c1 = list(P = diag(c(4,2,2,0)),r = c(2,-1,0,-1),b = 5),
       c2 = list(P = diag(c(2,2,2,2)),r = c(1,-1,1,-1),b = 8),
       c3 = list(P = diag(c(2,4,2,4)),r = c(-1,0,0,-1),b = 10))

# Solve the quadratic program using the primal-dual interior-point solver.
out <- ipsolver(x      = rep(0,4),
                obj    = function (x) quadf(x,H,u),
                grad   = function (x) list(g = c(H %*% x + u),H = H),
                constr = function (x) sapply(constraints,
                                        function (a) with(a,quadf(x,P,r,-b))),
                jac    = function (x, z) qpjacobian(x,z,constraints))
cat("Solution:\n")
options(digits = 4)
print(out$x)

