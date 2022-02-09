#' Wasserstein distribution-free statistic on S^1
#'
#' Computes Wasserstein statistic between two samples on S^1, parameterized as a periodic [0,1).
#'
#' @param x numerical vector of sample in [0,1).
#' @param y numerical vector of sample in [0,1).
#'
#' @return The statistic realization of samples x, y.
#'
#' @references [1] RABIN, J., DELON, J. and GOUSSEAU, Y. (2009). Transportation Distances on the Circle. Journal of Mathematical Imaging and Vision 41.
#' [2] RAMDAS, A., GARCIA, N. and CUTURI, M. (2015). On Wasserstein Two Sample Testing and Related Families of Nonparametric Tests. Entropy 19.
#' 
#' @examples
#' set.seed(10)
#' stat.s1(runif(50),runif(50))
#' 0.04916085
#' stat.s1(runif(50),as.numeric(circular::rvonmises(50,pi,1)/(2*pi)))
#' 0.1626842
#'  
#' @export

stat.s1<- function(x, y){
  
  fy <- ecdf(y)
  h <- function(t){(fy(quantile(x, t)) - t)}
  optimal_alpha <- pracma::integral(h,0,1)
  
  h_alpha <- function(t){(fy(quantile(x, t)) - t - optimal_alpha)^2}
  
  n <- length(x);  m <- length(y)
  return(n*m / (n + m) * pracma::integral(h_alpha, 0, 1))
  
}

