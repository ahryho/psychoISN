#' Title
#'
#' @param G 
#' @param meth 
#'
#' @return
#' @export
#'
#' @examples
euclidean_dist <- function(G, meth){
  nbNet <- length(G) #number of networks
  dist_mtrx <- matrix(ncol = nbNet, nrow = nbNet) # Distance matrix between graphs
  
  for(i in 1:nbNet){
    for(j in 1:nbNet){
      if (i == j)
        dist_mtrx[i,j] <- 0 else if(j > i){
        dist_mtrx[i,j] <- sqrt(sum((G[[i]] - G[[j]])^2))
      }
    }
  }
  
  dist_mtrx <- format(dist_mtrx, nbNet)
  
  return(dist_mtrx)
}

#' Title
#'
#' @param mat 
#' @param nbNet 
#'
#' @return
#' @export
#'
#' @examples
format <- function(mat, nbNet){
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  rownames(mat) <- seq(1, nbNet, 1)
  colnames(mat) <- seq(1, nbNet, 1)
  return(mat)
}