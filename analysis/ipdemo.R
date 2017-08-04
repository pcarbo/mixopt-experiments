# Demonstration of the primal-dual interior-point solver for fitting a
# logistic regression model. It computes estimates of the regression
# coefficients subject to an L1 penalization ("Lasso").

# SCRIPT PARAMETERS
# -----------------
n      <- 1000  # Number of data examples.
lambda <- 1/2   # L1 penalty strength.

# Ground-truth regression coefficients.
beta <- c(0,0,2,-4,0,0,-1,3) 

# LOAD PACKAGES AND FUNCTIONS
# ---------------------------
source("../code/ipsolver.R")

# GENERATE DATA SET
# -----------------
# This function returns the sigmoid function at x.
sigmoid <- function (x)
  1/(1 + exp(-x))

# Get the number of variables/features.
cat("Generating data set.\n")
p <- length(beta)
X <- matrix(rnorm(n*p),n,p)
y <- as.numeric(runif(n) < c(sigmoid(X %*% beta)))

# FIT MODEL
# ---------
# This function returns the logistic loss function with a penalty
# term given by the L1 norm of the regression coefficients. Here, it
# is assumed that the regression coefficients w are split into the
# positive and negative components so that all entries of w are
# positive. Argument a is the L1 penalty strength.
logisticl1.obj <- function (X, y, w, a) {
  u <- c(sigmoid(X %*% w))
  return(a*sum(w) - sum(y*log(u) + (1 - y)*log(1 - u)))
}

# This function returns the gradient and Hessian of the L1-penalized
# log-likelihood objective.
logisticl1.grad <- function (X, y, w, a) {
  u <- c(sigmoid(X %*% w))
  return(list(g = c(t(X) %*% (y - u)) + a,
              H = t(X) %*% diag(u*(1 - u)) %*% X))
}

# Minimize the penalized L1-penalized log-likelihood objective 
A   <- rbind(X,-X)
fit <- ipsolver(x      = rep(1,2*p),
                obj    = function (w) logisticl1.obj(A,y,w,lambda),
                grad   = function (w) logisticl1.grad(A,y,w,lambda),
                constr = function (w) -w,
                jac    = function (w) list(J = -diag(rep(1,p)),
                                           W = matrix(0,p,p)))
w  <- fit$w[1:p] - fit$w[-(1:p)]
