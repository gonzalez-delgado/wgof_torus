#' Two-sample goodness-of-fit on the circle
#'
#' Performs a two-sample goodness-of-fit test for measures supported on the circle, based on a
#' distribution-free Wasserstein statistic.
#'
#' @param x numerical vector of sample in [0,1).
#' @param y numerical vector of sample in [0,1).
#' @param sim_null The simulated null distribution of the statistic. It is introduced as 
#' an argument to avoid re-running the same simulation for each pair of equally-sized samples.
#' 
#' @return
#' \itemize{
#'   \item stat - The test statistic.
#'   \item pvalue - The test p-value.
#' }
#'
#' @references [1] RAMDAS, A., GARCIA, N. and CUTURI, M. (2015). On Wasserstein Two Sample Testing and Related Familites of Nonparametric Tests. Entropy 19.
#'
#'
#' @examples
#' 
#' n <- 30 # Sample size
#' 
#' # Simulate the statistic null distribution
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
#' x <- runif(n, 0, 1)
#' y <- runif(n, 0, 1)
#' twosample.test.s1(x, y, sim_free_null) 
#' 
#' x <- as.numeric(circular::rvonmises(n, pi, 1)/(2*pi))
#' y <- as.numeric(circular::rvonmises(n, pi, 0)/(2*pi))
#' twosample.test.s1(x, y, sim_free_null) 
#' 
#' @export

twosample.test.s1 <- function(x, y, sim_null){
  
  cut_point<-find.cutpoint.s1(x, y)$cut_point #Optimal origin
  
  x_left <- sort(x[which(x < cut_point)])
  x_right <- sort(x[which(x >= cut_point)])
  x_new <- c(x_right,x_left+1)
  x_new <- x_new - cut_point #Relocated measure on [0,1)
  
  y_left <- sort(y[which(y < cut_point)])
  y_right <- sort(y[which(y >= cut_point)])
  y_new <- c(y_right, y_left + 1)
  y_new <- y_new - cut_point #Relocated measure on [0,1)
  
  statistic <- stat.s1(x_new, y_new) #Ramdas statistic on [0,1)
  pv <- mean(sim_null > statistic) #Ramdas p-value
  
  return(list(stat = statistic, pvalue = pv))
}


