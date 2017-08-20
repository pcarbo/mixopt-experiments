# TO DO: Explain here what this script does.
#
# NOTES:
#
#   * This is a quadratic program with an equality constraint and an
#     inequality constraint.
#
#   * See problem #14 from H&S test problems.
#

# Load the packages and function definitions used in this demo.
library(Matrix)
source("../code/misc.R")
source("../code/ipsolver.R")

# DEFINE PROBLEM
# --------------
# TO DO.

# SOLVE PROBLEM USING IP METHOD
# -----------------------------
# Solve the quadratic program using the primal-dual interior-point
# solver.
out <- ipsolver(x      = rep(0,4),
                obj    = function (x) quadf(x,H,u),
                grad   = function (x) list(g = drop(H %*% x + u),H = H),
                constr = function (x) sapply(constraints,
                                        function (a) with(a,quadf(x,P,r,-b))),
                jac    = function (x, z) qpjacobian(x,z,constraints))
cat("Solution:\n")
print(out$x)
