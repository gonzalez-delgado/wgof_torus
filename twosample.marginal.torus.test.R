#' Marginal two-sample goodness-of-fit on the torus
#'
#' Performs a two-sample goodness-of-fit test for measures supported on the torus, by testing
#' the equality of their marginal laws on the circle and combining the results into a global p-value.
#'
#' @param sample_1 n x 2 matrix containing n observations in the two-dimensional flat torus, parameterized as the periodic [0,1) x [0,1).
#' @param sample_2 n x 2 matrix containing n observations in the two-dimensional flat torus, parameterized as the periodic [0,1) x [0,1).
#' @param sim_null The simulated null distribution of the circle test statistic. It NULL, the null distribution is simulated with the given parameters (very time consuming).
#' @param NR The number of replicas if simulation is required.
#' @param NC The number of cores if parallel simulation is required.
#' @param n The sample sizes of the simulated samples is simulation is required.
#' 
#' @return The p-value for the two-sample test on the torus.
#'
#' @examples
#' 
#' n <- 200 # Sample size
#' 
#' # Simulate the null distribution of the circle test statistic
#' sim_free_null <- sim.null.stat(500)
#' 
#' samp_1<-BAMBI::rvmcos(n, kappa1=1, kappa2=1, mu1=0, mu2=0)/(2*pi) # Bivariate von Mises distributions
#' samp_2<-BAMBI::rvmcos(n, kappa1=1, kappa2=1, mu1=0, mu2=0)/(2*pi)
#' twosample.marginal.torus(samp_1, samp_2, sim_free_null) 
#'
#' samp_1<-BAMBI::rvmcos(n ,kappa1=0, kappa2=0, mu1=0.5, mu2=0.5)/(2*pi)
#' samp_2<-BAMBI::rvmcos(n, kappa1=1, kappa2=1, mu1=0.5, mu2=0.5)/(2*pi)
#' twosample.marginal.torus(samp_1, samp_2, sim_free_null) 
#' 
#' @export

twosample.marginal.torus<-function(sample_1, sample_2, sim_null=NULL, NR = 500, NC = 1, n=30){
  
  p_1 <- twosample.test.s1(sample_1[,1], sample_2[,1], sim_null, NR, NC, n)$pvalue #First marginal p-value
  p_2 <- twosample.test.s1(sample_1[,2], sample_2[,2], sim_null, NR, NC, n)$pvalue #Second marginal p-value
  
  return(pbeta(min(p_1,p_2),1,2)) #p-value for the bi-dimensional test

}
