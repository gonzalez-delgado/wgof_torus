
#' Distance on the circle
#'
#' Distance between two points on the circle (periodic [0,1])
#'
#' @param x a number in [0,1].
#' @param y a number in [0,1].
#'
#' @return The distance on the circle between the point x and y
#'
#' @examples
#' set.seed(10)
#' dist.s1(runif(1),runif(1))
#' 0.2007097
#' dist.s1(runif(1),rbeta(1,1,2))
#' 0.245644
#' 
#' @export

dist.s1 <- function(x, y){
  
  dist <- sqrt(pmin(abs(x[1] - y[1]), 
                    1 - abs(x[1] - y[1]))^2)
  
  return(dist)
}
