# TO DO: Give overview of this analysis, and instructions for
# experimenting with different analysis parameters.

# DATA SIMULATION SETTINGS
# ------------------------
# Random number generator seed.
seed <- 1

# Number of data samples.
n <- 1e4

# The standard deviations and mixture weights used to simulate the data. 
sim <- list(s = c(0,   0.1, 0.2, 0.5),
            w = c(0.95,0.03,0.01,0.01))

# MODEL PARAMETERS
# ----------------
# The standard deviations of the normal mixture components. 
s <- c(0.01,10^(seq(-2,0,length.out = 39)))

# OPTIMIZATION SETTINGS
# ---------------------
# TO DO.

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
se <- rep(0.1,n)
x  <- datasim.norm(sim$w,sim$s,se)

# COMPUTE LIKELIHOOD MATRIX
# -------------------------
k <- length(s)
cat(sprintf("Computing the %d x %d conditional likelihood matrix.\n",n,k))
L <- condlikmatrix.norm(x,se,s)

# FIT MIXTURE MODEL
# -----------------
cat("Fitting model to data.\n")
out <- system.time(fit.em <- mixopt.em(L))
cat(sprintf("Model fitting took %0.2f seconds.\n",out["elapsed"]))
