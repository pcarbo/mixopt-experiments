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

# Return

# Quadratic objective function.
H <- diag(c(2,2,4,2))
u <- c(-5,-5,-21,7)

# Quadratic inequality constraints.
constraints <-
  list(c1 = list(P = diag(c(4,2,2,0)),r = c(2,-1,0,-1),b = 5),
       c2 = list(P = diag(c(2,2,2,2)),r = c(1,-1,1,-1),b = 8),
       c3 = list(P = diag(c(2,4,2,4)),r = c(-1,0,0,-1),b = 10))

# Solve the quadratic program using the primal-dual interior-point solver.
out <- ipsolver(x      = c(0,0,0,0),
                obj    = function (x) quadf(H,u,x)
                grad   = function (x) list(g = H %*% x + u,H = H),
                constr = function (x) sapply(constraints,
                                        function (a) with(a,quadf(x,P,r,-b))),
                jac    = qpjacobian(x,z))
cat("Solution:\n")
print(out$x)

    for i = 1:m
      J(i,:) = (P{i}*x + r{i})';
      W      = W + z(i)*P{i};
    end
