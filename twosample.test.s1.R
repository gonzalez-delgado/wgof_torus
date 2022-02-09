#' Two-sample goodness-of-fit on the circle
#'
#' Performs a two-sample goodness-of-fit test for measures supported on the circle, based on a
#' distribution-free Wasserstein statistic.
#'
#' @param x numerical vector of sample in [0,1).
#' @param y numerical vector of sample in [0,1).
#' @param sim_null The simulated null distribution of the statistic. If NULL, the distribution is simulated with the given parameters (very time consuming).
#' @param NR The number of replicas if simulation is required.
#' @param NC The number of cores if parallel simulation is required.
#' @param n The sample sizes of the simulated samples is simulation is required.
#' 
#' @return
#' \itemize{
#'   \item stat - The test statistic.
#'   \item pvalue - The test p-value.
#' }
#'
#' @references [1] RABIN, J., DELON, J. and GOUSSEAU, Y. (2009). Transportation Distances on the Circle. Journal of Mathematical Imaging and Vision 41.
#' [2] RAMDAS, A., GARCIA, N. and CUTURI, M. (2015). On Wasserstein Two Sample Testing and Related Families of Nonparametric Tests. Entropy 19.
#' 
#' @examples
#' 
#' n <- 50 # Sample size
#' 
#' # Simulate the statistic null distribution
#' NR <- 100
#' sim_free_null <- sim.null.stat(500,NC=3)
#' 
#' x <- runif(n, 0, 1)
#' y <- runif(n, 0, 1)
#' twosample.test.s1(x, y, sim_free_null) 
#' 
#' x <- as.numeric(circular::rvonmises(n, pi, 1)/(2*pi))
#' y <- as.numeric(circular::rvonmises(n, pi, 0)/(2*pi))
#' twosample.test.s1(x, y, sim_free_null) 
#' 
#' @export

twosample.test.s1 <- function(x, y, sim_null = NULL, NR = 500, NC = 1, n=30){
  
  if(is.null(sim_null)){
    
    cat('No null distribution given as an argument. Simulating with default parameters...\n')
    sim_null <- sim.null.stat(NR, NC, n)
    cat('Done.\n')
  }
  
  statistic <- stat.s1(x, y) #Ramdas statistic on [0,1)
  pv <- mean(sim_null > statistic) #Ramdas p-value
  
  return(list(stat = statistic, pvalue = pv))
}


