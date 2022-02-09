#' Wasserstein statistic null distribution
#'
#' Simulates Wasserstein statistic null distribution
#'
#' @param NR Number of replicas to simulate
#' @param NC Number of cores for parallel computation
#' @param n Sample size for the simulated samples
#'
#' @return A sample of size NR simulating the null distribution of the statistic.
#'
#' @references [1] RABIN, J., DELON, J. and GOUSSEAU, Y. (2009). Transportation Distances on the Circle. Journal of Mathematical Imaging and Vision 41.
#' [2] RAMDAS, A., GARCIA, N. and CUTURI, M. (2015). On Wasserstein Two Sample Testing and Related Families of Nonparametric Tests. Entropy 19.
#' 
#' @examples
#' sim.null.stat(100)
#'  
#' @export

sim.null.stat<- function(NR, NC = 1, n = 30){
  
  stat_rep <- function(i, m){
    
    x <- runif(m, 0, 1)
    y <- runif(m, 0, 1)
    fy <- ecdf(y)
    h <- function(t){(fy(quantile(x, t)) - t)}
    optimal_alpha <- pracma::integral(h,0,1)
    
    h_alpha <- function(t){(fy(quantile(x, t)) - t - optimal_alpha)^2}
    
    n <- length(x);  m <- length(y)
    return(n*m / (n + m) * pracma::integral(h_alpha, 0, 1))
  
    }
  
  clus <- parallel::makeCluster(NC)
  doParallel::registerDoParallel(clus)
  
  sim <- parallel::parApply(clus, X = as.matrix(1:NR), MARGIN = 1, FUN = stat_rep, m=n)
  parallel::stopCluster(clus)
  
  return(sim)
}

