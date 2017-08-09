# TO DO: Give overview of this analysis, and instructions for
# experimenting with different analysis parameters.

# DATA SIMULATION SETTINGS
# ------------------------
# Random number generator seed.
seed <- 1

# Number of data samples.
n <- 200

# The standard deviations and mixture weights used to simulate the data. 
sim <- list(s = c(0,   0.1, 0.2, 0.5),
            w = c(0.95,0.03,0.01,0.01))

# MODEL PARAMETERS
# ----------------
# The standard deviations of the normal mixture components. 
s <- c(0.01,10^(seq(-2,0,length.out = 39)))

# LOAD PACKAGES AND FUNCTIONS
# ---------------------------
# Load the data simulation functions, likelihood computation functions
# and optimization algorithms into the R environment.
library(Matrix)
library(scales)
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
out <- system.time(fit.em <- mixopt.em(L,tol = 1e-4))
cat(sprintf("Model fitting took %0.2f seconds.\n",out["elapsed"]))

# FIT MIXTURE MODEL USING IP METHOD
# ---------------------------------
# Compare against
#
#   fit.ip <- REBayes:KWDual(L,rep(1,k),rep(1,n))
#
cat("Fitting model using interior-point algorithm.\n")
out <- system.time(fit.ip <- mixopt.dualip(L))
cat(sprintf("Model fitting took %0.2f seconds.\n",out["elapsed"]))

# PLOT OPTIMIZATION RESULTS
# -------------------------
adjust.plot <- function (p)
  p + theme_cowplot(font_size = 9) +
    theme(plot.title = element_text(face = "plain"))

# Show the (maximum) change in the mixture weights vs. time running
# the EM algorithm.
m  <- length(fit.em$maxd)
i  <- 2:(m-1)
p1 <- ggplot(data.frame(time = fit.em$timing[i,"elapsed"],
                        maxd = fit.em$maxd[i]),
             aes(x = time,y = maxd)) +
  geom_line(col = "darkorange",size = 0.5) +
  geom_point(col = "darkorange",shape = 20) +
  scale_y_continuous(breaks = 10^seq(-4,-1),trans = "log10") +
  labs(x     = "elapsed time (seconds)",
       y     = "max. change in solution",
       title = "EM algorithm")

# Show the (maximum) change in the mixture weights vs. time running
# the interior-point method.
m  <- length(fit.ip$maxd)
i  <- 2:(m-1)
p2 <- ggplot(data.frame(time = fit.ip$timing[i,"elapsed"],
                        maxd = fit.ip$maxd[i]),
             aes(x = time,y = maxd)) +
  geom_line(col = "darkblue",size = 0.5) +
  geom_point(col = "darkblue",shape = 20) +
  scale_y_continuous(breaks = 10^(-3:1),trans = "log10") +
  labs(x     = "elapsed time (seconds)",
       y     = "max. change in solution",
       title = "IP algorithm")

# Show the value of the (primal) objective function vs. elapsed time
# running the EM algorithm.
m  <- length(fit.em$obj)
i  <- 2:(m-1)
p3 <- ggplot(data.frame(time = fit.em$timing[i,"elapsed"],
                        y    = fit.em$obj[i] - fit.ip$obj),
                        aes(x = time,y = y)) +
  geom_line(col = "darkorange",size = 0.5) +
  geom_point(col = "darkorange",shape = 20) +
  scale_y_continuous(breaks = 10^(-1:1),trans = "log10") +
  labs(x     = "elapsed time (seconds)",
       y     = "distance from primal min.",
       title = "EM algorithm")

# Show the value of the (dual) objective function vs. elapsed time
# running the interior-point method.
m  <- length(fit.ip$ipsolver$obj)
i  <- 2:(m-1)
p4 <- ggplot(data.frame(time = fit.ip$timing[i,"elapsed"],
                        y = fit.ip$ipsolver$obj[i] - min(fit.ip$ipsolver$obj)),
                        aes(x = time,y = y)) +
  geom_line(col = "darkorange",size = 0.5) +
  geom_point(col = "darkorange",shape = 20) +
  scale_y_continuous(breaks = 10^seq(-8,2,2),trans = "log10") +
  labs(x     = "elapsed time (seconds)",
       y     = "distance from dual min.",
       title = "IP algorithm")

# Draw all four plots.
print(plot_grid(adjust.plot(p1),adjust.plot(p3),
                adjust.plot(p2),adjust.plot(p4)))

