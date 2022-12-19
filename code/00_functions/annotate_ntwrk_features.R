#' Title
#'
#' @param net 
#' @param anno_gr 
#' @param out_fn 
#'
#' @return
#' @export
#'
#' @examples
annotate_ntwk_cpgs <- function(net, anno_gr, out_fn = NULL){
  net_features     <- rownames(net)[rownames(net) %like% "cg"]
  anno_net_feat_gr <- anno_gr[names(anno_gr) %in% net_features]
  
  if(!is.null(out_fn))
    saveRDS(anno_net_feat_gr, out_fn)
  
  return(anno_net_feat_gr)
}

#' Title
#'
#' @param net 
#' @param anno_gr 
#' @param out_fn 
#'
#' @return
#' @export
#'
#' @examples
annotate_ntwk_snps <- function(net, anno_gr, out_fn = NULL){
  net_features     <- rownames(net)[rownames(net) %like% "rs"]
  anno_net_feat_gr <- anno_gr[names(anno_gr) %in% net_features]
  
  if(!is.null(out_fn))
    saveRDS(anno_net_feat_gr, out_fn)
  
  return(anno_net_feat_gr)
}

#' Title
#'
#' @param feature_gr 
#' @param chromhmm_states_gr 
#' @param out_fn 
#'
#' @return
#' @export
#'
#' @examples
annotate_chromHMM <- function(feature_gr, chromhmm_states_gr, out_fn){
  
  annotated_gr <- annotate_regions(
    regions = GRanges(seqnames = seqnames(feature_gr), ranges = ranges(feature_gr)),
    annotations = chromhmm_states_gr,
    ignore.strand = TRUE,
    quiet = FALSE)
  
  annotated_df         <- data.frame(annotated_gr) %>% setDT()
  annotated_df[["ID"]] <- names(annotated_gr)
  annotated_df         <- annotated_df[, .(ID, seqnames, start, annot.type, annot.code)] %>% 
    unique() %>% select(ID, chr = seqnames, pos = start, annot.type, annot.code)
  
  fwrite(annotated_df,
         out_fn,
         quote = F, row.names = F, sep = "\t")
  
  return(annotated_df)
}
