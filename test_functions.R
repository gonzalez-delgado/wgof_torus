library(pracma)
library(transport)
library(dplyr)
library(som.nn)

# 1D -marginal test based on Wasserstein distance

# Required functions to perform the test:

# Computes distance on the torus (periodic [0,1]x[0,1])
dist.torus.mod <- function(coor){
  dist <- sqrt(pmin(abs(coor[1] - coor[3]), 
                    1 - abs(coor[1] - coor[3]))^2 + 
                 pmin(abs(coor[2] - coor[4]), 
                      1-abs(coor[2] - coor[4]))^2)
  return(dist)
}

#Computes Ramdas statistic
stat_1d <- function(x,y){
  fy <- ecdf(y)
  h <- function(t){
    (fy(quantile(x,t))-t)^2
  }
  n <- length(x)
  m <- length(y)
  
  return(n*m / (n + m) * pracma::integral(h,0,1))
}

# Computes Ramdas p-value by simulating the statistic distribution under the null
test_1d <- function(x, y, NR = 1000, NS = 1000){
  # Simulate test statistic distribution under the null (P=Q), which is free (Ramadas et al.)
  sim_dist <- replicate(NR, {
    x1 <- runif(NS, 0, 1)
    x2 <- runif(NS, 0, 1)
    stat_1d(x1, x2)
  })
  
  statistic <- stat_1d(x, y)
  pv <- mean(sim_dist > statistic)
  
  return(list(stat = statistic, pvalue = pv))
}

#Finds the cut-point on the circle and sends the problem to the real line
transport_circ <- function(x, y, by_seq = 0.01){
  
  n <- min(length(x), length(y))
  if(n == length(x)){
    y <- sample(y,n)
  } else{
    x <- sample(x,n)
  } #Same sample size
  
  coord <- expand.grid(1:n, 1:n)
  colnames(coord) <- c('from', 'to')
  coord$x <- x[coord$from]
  coord$xf <- 0
  coord$y <- y[coord$to]
  coord$yf <- 0
  
  coord$cost <- as.vector(apply(X = coord[,3:6], FUN = dist.torus.mod, MARGIN = 1))
  mat <- matrix(coord$cost, nrow = n, ncol = n) #Cost matrix
  
  a <- rep(1/n, n)
  b <- rep(1/n, n)
  ot <- transport::transport(a, b, p = 2, costm = mat, method = "networkflow") # Optimal transport map
  
  lj <- left_join(coord, ot, by = "from")
  ot_map <- lj %>% filter(to.x == to.y)
  
  cutpt <- seq(0, 1, by = by_seq)
  cutpt <- sample(cutpt,size = length(cutpt))
  cat('Looking for a cutpoint\n')
  
  for (i in 1:length(cutpt)) {
    coor <- as.data.frame(cbind(
      matrix(rep(c(cutpt[i], 0), n), nrow = n, byrow = TRUE), 
      matrix(c(ot_map$x, rep(0, n)), nrow = n, byrow = FALSE),
      matrix(c(ot_map$y, rep(0, n)), nrow = n, byrow = FALSE)))
    colnames(coor) <- c('seqpt', 'V2', 'x', 'V2x', 'y', 'V2y')
    coor$dis_x <- apply(coor[, 1:4], 1, dist.torus.mod)
    coor$dis_y <- apply(coor[, c(1:2,5:6)], 1, dist.torus.mod)
    coor$opt_dis <- ot_map$cost
    coor$dis_seqpt <- coor$dis_x + coor$dis_y
    coor$cut <- coor$dis_seqpt == coor$opt_dis
    if (sum(coor$cut) == 0) {
      cat('Cut point has been found.\n')
      i_opt <- i
      return(list(opt_plan = ot, cut_point = cutpt[i_opt], x = x, y = y))
    } else {
      if (i == length(cutpt)) {
        stop('No cut point has been found. Candidate set must be thinned.')
      }
      next
    }
  }
}

# Performs the 1D-marginal test using the previous functions
test_1d_circ <- function(x, y, by_seq = 0.01, NR = 500, NS = 1000){
  
  tc <- transport_circ(x, y, by_seq)
  
  x_left <- sort(tc$x[which(tc$x < tc$cut_point)])
  x_right <- sort(tc$x[which(tc$x >= tc$cut_point)])
  x_new <- c(x_right,x_left+1)
  
  y_left <- sort(tc$y[which(tc$y < tc$cut_point)])
  y_right <- sort(tc$y[which(tc$y >= tc$cut_point)])
  y_new <- c(y_right, y_left + 1)
  
  return(test_1d(x_new, y_new, NR, NS))
}

################################################################################

# Upper bound test

# Compute empirical Wasserstein distance on the torus
wd_torus <- function(x1, y1, x2, y2, p = 2, ot_method = "networkflow") {
  
  grid <- expand.grid(1:length(x1), 1:length(x2))
  grid$x1 <- x1[grid[, 1]]
  grid$y1 <- y1[grid[, 1]]
  grid$x2 <- x2[grid[, 2]]
  grid$y2 <- y2[grid[, 2]]
  
  grid$dist <- as.vector(apply(grid[, 3:6], FUN = dist.torus.mod, MARGIN = 1))
  mat <- matrix(grid$dist, nrow = length(x1), ncol = length(x2)) #Cost matrix
  
  n <- length(x1)
  m <- length(x2)
  dist <- transport::wasserstein(rep(1,n), rep(1,m), 
                                 p, costm = mat, method = ot_method) 
  return(dist^p) # Statistic value
  }

# Upper bound
u_bound <- function(t, n, m){
  exp(-8 * t^2 * m * n / (n+m))
}
