#' Optimal origin on the circle
#'
#' Finds an optimal origin for equally-sized finite samples on the circle
#'
#' @param x numerical vector of sample in [0,1).
#' @param y numerical vector of sample in [0,1).
#' @param message (logical) whether to inform of the computation progress.
#' 
#' @return
#' \itemize{
#'   \item opt_plan - The optimal transport plan between both empirical measures on the circle.
#'   \item cut_point - An admissible optimal origin to relocate both samples on the real interval [0,1).
#' }
#'
#' @references [1] RABIN, J., DELON, J. and GOUSSEAU, Y. (2009). Transportation Distances on the Circle. Journal of Mathematical Imaging and Vision 41.
#' 
#' @examples
#' find.cutpoint.s1(runif(50),runif(50))
#' find.cutpoint.s1(runif(50),rbeta(50,1,2))
#'  
#' @export

find.cutpoint.s1<- function(x, y, messages = TRUE){
  
  n <- min(length(x), length(y))
  if(n == length(x)){
    y <- sample(y,n)
  } else{
    x <- sample(x,n)
  } #Force same sample size
  
  coord <- expand.grid(1:n, 1:n)
  colnames(coord) <- c('from', 'to')
  coord$x <- x[coord$from]
  coord$y <- y[coord$to]
  
  if(messages){cat('Solving OT problem...\n')}
  mat<-proxy::dist(x=x, y=y, method=dist.s1, diag=TRUE) #Cost matrix
  a <- rep(1/n, n)
  b <- rep(1/n, n)
  ot <- transport::transport(a, b, p = 2, costm = mat, method = "networkflow") # Optimal transport map
  
  lj <- dplyr::left_join(coord, ot, by = "from") 
  ot_map <- dplyr::filter(lj,to.x == to.y) 
  if(messages){cat('Looking for a cutpoint...\n')}
  
  i_list<-list(); k<-1
  
  for(i in 1:nrow(ot_map)){
    
    left <- min(ot_map$x[i], ot_map$y[i])
    right <- max(ot_map$x[i], ot_map$y[i])
    dis_i1 <- abs(left - right)
    dis_i2 <- 1 - dis_i1
    
    if(min(dis_i1, dis_i2) == dis_i1){
      
      i_list[[k]] <- sets::interval(l = left, r = right)}else{
        
        i_list[[k]] <- sets::interval(0, left, bounds = "[)")
        k<-k+1
        i_list[[k]] <- sets::interval(right, 1, bounds = "(]")
        
      }
    k<-k+1
  }
  
  i_union <- sets::interval_union(i_list)
  #The sub-interval of [0,1) containing all the admissible optimal origins
  admisible_cutpoints <- sets::interval_intersection(sets::interval(0,1), sets::interval_complement(i_union))
  if(messages){cat('Done.\n')}
  
  #Choose an optimal origin among all the admissible
  admisible_cutpoints <- stringr::str_split(as.character(admisible_cutpoints), pattern = ' U')[[1]][1]
  adm_interval <- as.numeric(as.vector(stringr::str_split(substring(admisible_cutpoints, 2, nchar(admisible_cutpoints)-1), pattern = ',')[[1]]))
  cutpoint <- runif(1, adm_interval[1], adm_interval[2]) #An admissible optimal origin
  
  return(list(opt_plan = ot, cut_point = cutpoint))
}
