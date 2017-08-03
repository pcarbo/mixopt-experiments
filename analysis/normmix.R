# TO DO: Give overview of this analysis, and instructions for
# experimenting with different analysis parameters.

# ANALYSIS PARAMETERS
# -------------------
seed <- 1     # Random number generator seed.
n    <- 5000  # Number of data samples.

# The standard deviations and mixture weights used to simulate the
# data. 
s <- c(0,   0.1, 0.2, 0.5)
w <- c(0.95,0.03,0.01,0.01)

# LOAD PACKAGES AND FUNCTIONS
# ---------------------------
# Load the data simulation functions, likelihood computation functions
# and optimization algorithms into the R environment.
source("../code/datasim.R")
source("../code/likelihood.R")
source("../code/mixopt.R")

# GENERATE DATA SET
# -----------------
# Set the random number generator seed.
set.seed(1)

# Simulate a data set with n samples, in which the standard errors are ...
cat(sprintf("Simulating data set with %d observations.\n",n))
k  <- length(w)
se <- rep(0.1,n)
x  <- datasim.norm(w,s,se)

# COMPUTE LIKELIHOOD MATRIX
# -------------------------
cat(sprintf("Computing the %d x %d conditional likelihood matrix.\n",
            n,length(w)))
L <- condlikmatrix.norm(x,se,s)

# FIT MIXTURE MODEL
# -----------------
cat("Fitting model to data.\n")
fit.em <- mixopt.em(L)

