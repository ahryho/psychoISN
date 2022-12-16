#' Title
#'
#' @param mtrx 
#'
#' @return
#' @export
#'
#' @examples
flatten_matrix <- function(mtrx) {
  if (is.null(mtrx)) return (NULL)
  
  ut_mtrx <- upper.tri(mtrx)
  
  flat_mtrx <- data.frame(
      row    = rownames(mtrx)[row(mtrx)[ut_mtrx]],
      column = rownames(mtrx)[col(mtrx)[ut_mtrx]],
      weight = (mtrx)[ut_mtrx])
  
  return(flat_mtrx)
}
