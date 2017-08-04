# Demonstration of the primal-dual interior-point solver for fitting a
# logistic regression model. It computes estimates of the regression
# coefficients subject to an L1 penalization ("Lasso").

# SCRIPT PARAMETERS
# -----------------
n      <- 1000  # Number of data examples.
lambda <- 1/2   # L1 penalty strength.

# Ground-truth regression coefficients.
beta <- c(0,0,2,-4,0,0,-1,3) 

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
# This function computes the logistic loss function with a penalty
# term given by the L1 norm of the regression coefficients. Here, it
# is assumed that the regression coefficients w are split into the
# positive and negative components so that all entries of w are
# positive. Argument a is the L1 penalty strength.
logisticl1.obj <- function (X, y, w, a) {
  u <- c(sigmoid(X %*% w))
  return(a*sum(w) - sum(y*log(u) + (1 - y)*log(1 - u)))
}

# This function computes 
logisticl1.grad <- function (X, y, w, a) {
  u <- c(sigmoid(X %*% w))
  return(g = c(t(X) %*% (y - u)) + a,
         H = t(X) %*% diag(u*(1 - u)) %*% X)
}

# Estimate the regression coeffiicents subject to an L1 penalty.
## P      = [A -A];
## x      = ones(2*m,1);  # The inital point.
## z      = ones(2*m,1);  # A dummy variable; does nothing.
## data   = { P y lambda };
## x      = ipsolver(x,@(x)logisticl1(x,z,data,'objective'),...
## 		  @(x)logisticl1(x,z,data,'gradient'),...
## 		  @(x)logisticl1(x,z,data,'constraints'),...
## 		  @(x,z)logisticl1(x,z,data,'jacobian'),...
## 		  'steepest',1e-4,100,true);
## w      = x(1:m) - x(m+1:end);
## fprintf('\nSolution:\n');
## disp(w);
