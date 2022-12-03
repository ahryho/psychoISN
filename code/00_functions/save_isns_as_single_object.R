#' Title
#'
#' @param fn_lst 
#' @param treatment 
#'
#' @return
#' @export
#'
#' @examples
save_isns_as_single_object <- function(fn_lst, treatment, ofile = NULL){
  isns_lst <- lapply(fn_lst, function(idx){
    isn_obj    <- readRDS(idx)
    isn        <- isn_obj$network
    membership <- data.frame(network_id = isn_obj$DNA_ID, treatment = treatment) 
    return (list(isn, membership))
  }) 
  
  isns_nets        <- lapply(isns_lst, function(isn) isn[1])
  isns_memberships <- lapply(isns_lst, function(isn) isn[2]) %>% bind_rows()
  isns_lst         <- list(isns_nets, isns_memberships)
  
  if(!is.null(ofile)) saveRDS(isns_lst, ofile)
  
  return(dex_isns_lst)
}
