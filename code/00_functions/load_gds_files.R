
#' Title
#'
#' @param gds_fn 
#' @param is_ld 
#'
#' @return
#' @export
#'
#' @examples
LoadGenotype <- function(gds_fn, is_ld = F){
  snps_mtrx = tryCatch(
    {
      snps_gds  <- snpgdsOpen(gds_fn)
      
      # snpset  <- snpgdsLDpruning(snps_gds, ld.threshold=0.2, maf = 0.05, slide.max.bp = 1000000L, slide.max.n = 100, method = "corr")
      
      geno_obj  <- snpgdsGetGeno(snps_gds, with.id = T)
      snp_chr   <- read.gdsn(index.gdsn(snps_gds, "snp.chromosome"))
      
      snpgdsClose(snps_gds)
      
      snp_mtrx             <- geno_obj$genotype
      colnames(snp_mtrx)   <- geno_obj$snp.id
      rownames(snp_mtrx)   <- geno_obj$sample.id
      
      return(snp_mtrx)
    },
    error = function(cond){
      # message(cond)
      
      snps_gds_file         <- openfn.gds(gds_fn)
      
      snps_mtrx             <- read.gdsn(index.gdsn(snps_gds_file, "genotype"))
      colnames(snps_mtrx)   <- read.gdsn(index.gdsn(snps_gds_file, "snp_id"))
      rownames(snps_mtrx)   <- read.gdsn(index.gdsn(snps_gds_file, "sample_id"))
      
      closefn.gds(snps_gds_file)
      
      return(snps_mtrx)
      }
  )
  return(snps_mtrx)
}

#' Title
#'
#' @param gds_fn 
#' @param is_mad 
#'
#' @return
#' @export
#'
#' @examples
LoadMethyl <- function(gds_fn, is_mad = F){
  dnam_gds_file         <- openfn.gds(gds_fn)
  
  dnam_mtrx             <- read.gdsn(index.gdsn(dnam_gds_file, "beta_mtrx"))
  colnames(dnam_mtrx)   <- read.gdsn(index.gdsn(dnam_gds_file, "cpg_id"))
  rownames(dnam_mtrx)   <- read.gdsn(index.gdsn(dnam_gds_file, "sample_id"))
  
  closefn.gds(dnam_gds_file)
  
  return(dnam_mtrx)
}
