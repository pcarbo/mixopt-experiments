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
library(ggplot2)
library(cowplot)
source("../code/datasim.R")
source("../code/likelihood.R")
source("../code/ipsolver.R")
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

# FIT MIXTURE MODEL USING EM ALGORITHM
# ------------------------------------
cat("Fitting model using EM.\n")
out <- system.time(fit.em <- mixopt.em(L))
cat(sprintf("Model fitting took %0.2f seconds.\n",out["elapsed"]))

# FIT MIXTURE MODEL USING IP METHOD
# ---------------------------------
cat("Fitting model using interior-point algorithm.\n")
fit.ip <- 
                                        

# PLOT OPTIMIZATION RESULTS
# -------------------------
# Show the maximum change in the mixture weights at each iteration of
# the EM algorithm. Here, average r is the amount of time elapsed per
# iteration.
m  <- length(fit.em$maxd)
i  <- 2:(m-1)
r  <- out["elapsed"] / m
p1 <- ggplot(data.frame(time = r*i,maxd = fit.em$maxd[i]),
             aes(x = time,y = maxd)) +
    geom_line(col = "darkorange",size = 1) + 
    theme_cowplot(font_size = 10) +
    scale_y_log10() +
    labs(x     = "elapsed time (seconds)",
         y     = "max. change in solution")

# Show the value of the objective function at each iteration of the EM
# algorithm.
p2 <- ggplot(data.frame(time = r*i,y = fit.em$obj[i] - min(fit.em$obj)),
                        aes(x = time,y = y)) +
    geom_line(col = "darkorange",size = 1) +
    scale_y_log10() +
    theme_cowplot(font_size = 10) +
    labs(x     = "elapsed time (seconds)",
         y     = "distance from minimum")

# Draw the two plots side-by-side.
plot_grid(p1,p2)

