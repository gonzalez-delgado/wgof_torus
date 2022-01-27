#' Wasserstein distribution-free statistic on [0,1)
#'
#' Computes Ramdas Wasserstein statistic between two samples on [0,1)
#'
#' @param x numerical vector of sample in [0,1).
#' @param y numerical vector of sample in [0,1).
#'
#' @return The statistic realization of samples x, y.
#'
#' @references [1] RAMDAS, A., GARCIA, N. and CUTURI, M. (2015). On Wasserstein Two Sample Testing and Related Familites of Nonparametric Tests. Entropy 19.
#'
#' @examples
#' set.seed(10)
#' stat.s1(runif(50),runif(50))
#' 0.07188233
#' stat.s1(runif(50),rbeta(50,1,2))
#' 1.471972
#'  
#' @export

stat.s1<- function(x, y){
  
  fy <- ecdf(y)
  h <- function(t){
    (fy(quantile(x,t)) - t)^2
  }
  
  n <- length(x);  m <- length(y)
  return(n*m / (n + m) * pracma::integral(h, 0, 1))
  
}

