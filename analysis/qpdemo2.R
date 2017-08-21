# TO DO: Explain here what this script does.
#
# NOTES:
#
#   * This is a quadratic program with an equality constraint and an
#     inequality constraint.
#
#   * See problem #14 (QQR-T1-4) from H&S test problems.
#
#   * The minimium occurs at x1 = (sqrt(7) - 1)/2, x2 = (sqrt(7) + 1)/4.
#
#   * Starting point must be (primal) feasible.
#

# Load the packages and function definitions used in this demo.
library(Matrix)
source("../code/misc.R")
source("../code/ipsolver.R")

# DEFINE PROBLEM
# --------------
# This function returns the objective function at x.
obj <- function (x)
  (x[1] - 2)^2 + (x[2] - 1)^2

# This function returns the gradient and Hessian of the objective
# function at x.
grad <- function (x)
  list(g = 2*c(x[1] - 2,x[2] - 1),
       H = diag(c(2,2)))

# This function returns the inequality constraint function.
constr <- function (x)
  x[1]^2/4 + x[2]^2 - 1

# This function computes the m x n Jacobian of the inequality
# constraint function, and the n x n Hessian of the Lagrangian (minus
# the Hessian of the objective).
jac <- function (x, z)
  list(J = t(matrix(c(0.5*x[1],2*x[2]))),
       W = z*diag(0.5,2))
    
# SOLVE PROBLEM USING IP METHOD
# -----------------------------
# Solve the quadratic program using the primal-dual interior-point
# solver.
out <- ipsolver(x = c(0.5,0.75),A = t(matrix(c(1,-2))),b = -1,
                obj = obj,grad = grad,constr = constr,jac = jac)
cat("Solution:\n")
print(out$x)
