library(uniformly)
library(MASS)

# source("test_functions.R")
################################################################################
# Reproducible examples

# Set sample size
ns <- 500

##################################
# H0: Two uniform laws on the torus
samp_1 <- uniformly::runif_in_cube(ns, 2, 0.5, 0.5)
samp_2 <- uniformly::runif_in_cube(ns, 2, 0.5, 0.5)

# 1D-marginal test
test_x <- test_1d_circ(samp_1[, 1], samp_2[, 1], by_seq = 0.001)
test_y <- test_1d_circ(samp_1[, 2], samp_2[, 2], by_seq = 0.001)
pv_1d <- min(1, 2*min(test_x$pvalue, test_y$pvalue)) # 1D-test p-value

# Upper bound
stat <- wd_torus(samp_1[, 1],
                 samp_1[, 2],
                 samp_2[, 1],
                 samp_2[, 2],
                 p = 2)
pv_ub <- u_bound(stat, ns, ns)

#Both p-values
cat(paste0('\t 1D-marginals\t Upper bound\n p-value \t', pv_1d, '\t', pv_ub, '\n'))

##################################
# H1: Uniform and gaussian on the torus
samp_1 <- runif_in_cube(ns, 2, 0.5, 0.5)
myCov <- matrix(c(0.01, 0.0001, 0.0001, 0.01), nrow = 2, byrow = TRUE)
samp_2 <- mvrnorm(n = ns, mu = c(0.5, 0.5), Sigma = myCov)
samp_2 <- samp_2[which(samp_2[, 1] < 1 & 
                           samp_2[, 1] > 0 & 
                           samp_2[, 2] < 1 & 
                           samp_2[, 2] > 0),] # Limit to (0,1)x(0,1)
ind <- sample(1:nrow(samp_2), size = ns - nrow(samp_2)) # Variance is set in order to keep points inside (0,1)x(0,1). If some points fall outside those limits, they are replaced by one of the inside points.
samp_2 <- rbind(samp_2, samp_2[ind, ])

#1D-marginal test
test_x <- test_1d_circ(samp_1[, 1], samp_2[, 1], by_seq = 0.001)
test_y <- test_1d_circ(samp_1[, 2], samp_2[, 2], by_seq = 0.001)
pv_1d <- min(1, 2*min(test_x$pvalue, test_y$pvalue)) #1D-test p-value

# Upper bound
stat <- wd_torus(samp_1[, 1], samp_1[, 2], samp_2[, 1], samp_2[, 2], p = 2)
pv_ub <- u_bound(stat, ns, ns)

# Both p-values
cat(paste0('\t 1D-marginals\t Upper bound\n p-value \t', pv_1d, '\t', pv_ub, '\n'))

##################################
#H4: Two different gaussians with equal marginals
myCov <- matrix(c(0.02, 0.019, 0.019, 0.02), nrow = 2, byrow = TRUE)
samp_1 <- mvrnorm(n = ns, mu = c(0.5, 0.5), Sigma = myCov)
samp_1 <- samp_1[which(samp_1[, 1] < 1 & 
                           samp_1[, 1] > 0 & 
                           samp_1[, 2] < 1 & 
                           samp_1[, 2] > 0), ] # Limit to (0,1)x(0,1)
ind <- sample(1:nrow(samp_1), size = ns - nrow(samp_1)) # Variance is set in order to keep points inside (0,1)x(0,1). If some points fall outside those limits, they are replaced by one of the inside points.
samp_1 <- rbind(samp_1, samp_1[ind, ]) 

myCov <- matrix(c(0.02, -0.019, -0.019, 0.02), nrow = 2, byrow = TRUE)
samp_2 <- mvrnorm(n = ns, mu = c(0.5, 0.5), Sigma = myCov)
samp_2 <- samp_2[which(samp_2[, 1] < 1 & 
                           samp_2[, 1] > 0 & 
                           samp_2[, 2] < 1 & 
                           samp_2[, 2] > 0), ] #Limit to (0,1)x(0,1)
ind <- sample(1:nrow(samp_2), size = ns-nrow(samp_2)) # Variance is set in order to keep points inside (0,1)x(0,1). If some points fall outside those limits, they are replaced by one of the inside points.
samp_2 <- rbind(samp_2, samp_2[ind, ])

# 1D-marginal test
test_x <- test_1d_circ(samp_1[, 1], samp_2[, 1], by_seq = 0.001)
test_y <- test_1d_circ(samp_1[, 2], samp_2[, 2], by_seq = 0.001)
pv_1d <- min(1, 2*min(test_x$pvalue, test_y$pvalue)) # 1D-test p-value

# Upper bound
stat <- wd_torus(samp_1[, 1], samp_1[, 2], samp_2[, 1], samp_2[, 2], p = 2)
pv_ub <- u_bound(stat, ns, ns)

# Both p-values
cat(paste0('\t 1D-marginals\t Upper bound\n p-value \t', pv_1d, '\t', pv_ub, '\n'))
