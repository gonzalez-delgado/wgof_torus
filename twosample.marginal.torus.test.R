#' Marginal two-sample goodness-of-fit on the torus
#'
#' Performs a two-sample goodness-of-fit test for measures supported on the torus, by testing
#' the equality of their marginal laws on the circle and combining the results into a global p-value.
#'
#' @param sample_1 n x 2 matrix containing n observations in the two-dimensional flat torus, parametrized as the periodic [0,1) x [0,1).
#' @param sample_2 n x 2 matrix containing n observations in the two-dimensional flat torus, parametrized as the periodic [0,1) x [0,1).
#' @param sim_null The simulated null distribution of the circle test statistic. It is introduced as 
#' an argument to avoid re-running the same simulation for each pair of equally-sized samples.
#' 
#' @return The p-value for the two-sample test on the torus.
#'
#' @examples
#' 
#' n <- 30 # Sample size
#' 
#' # Simulate the null distribution of the circle test statistic
#' NR <- 100 
#' sim_free_null <- pbapply::pbreplicate(NR, {
#'  
#'  x <- runif(n, 0,1)
#'  y <- runif(n, 0,1)
#'  
#'  cut_point<-find.cutpoint.s1(x, y, messages = FALSE)$cut_point #Optimal origin
#'  
#'  x_left <- sort(x[which(x < cut_point)])
#'  x_right <- sort(x[which(x >= cut_point)])
#'  x_new <- c(x_right,x_left+1)
#'  x_new <- x_new - cut_point 
#'  
#'  y_left <- sort(y[which(y < cut_point)])
#'  y_right <- sort(y[which(y >= cut_point)])
#'  y_new <- c(y_right, y_left + 1)
#'  y_new <- y_new - cut_point 
#'  
#'  statistic <- stat.s1(x_new, y_new) 
#'  
#'  statistic
#' })
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

twosample.marginal.torus<-function(sample_1, sample_2, sim_null){
  
  p_1 <- twosample.test.s1(sample_1[,1], sample_2[,1], sim_null)$pvalue #First marginal p-value
  p_2 <- twosample.test.s1(sample_1[,2], sample_2[,2], sim_null)$pvalue #Second marginal p-value
  
  return(pbeta(min(p_1,p_2),1,2)) #p-value for the bi-dimensional test

}