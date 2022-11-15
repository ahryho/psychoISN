
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
      print("Loading Genotype  matrix... ", quote = F)
      
      snps_gds  <- snpgdsOpen(gds_fn)
      
      # snpset  <- snpgdsLDpruning(snps_gds, ld.threshold=0.2, maf = 0.05, slide.max.bp = 1000000L, slide.max.n = 100, method = "corr")
      
      geno_obj  <- snpgdsGetGeno(snps_gds, with.id = T)
      snp_chr   <- read.gdsn(index.gdsn(snps_gds, "snp.chromosome"))
      
      snpgdsClose(snps_gds)
      
      snps_mtrx             <- geno_obj$genotype
      colnames(snps_mtrx)   <- geno_obj$snp.id
      rownames(snsp_mtrx)   <- geno_obj$sample.id
      
      print("Genotype  matrix has been loaded", quote = F)
      print(paste0("DNAm matrix: ", nrow(snps_mtrx), " x ", ncol(snps_mtrx)), quote = F)
      return(snps_mtrx)
    },
    error = function(cond){
      # message(cond)
      
      snps_mtrx = tryCatch(
        {
          snps_gds_file         <- openfn.gds(gds_fn)
          
          snps_mtrx             <- read.gdsn(index.gdsn(snps_gds_file, "genotype"))
          
          # If node "snp_id" exists, load it, otherwise, load "snp.id"
          # The "snp_id" node name is relevant to data splitted by chrom, 
          # meanwhile "snp.id" node name is relevent for full unsplitted genome 
          colnames(snps_mtrx)   <- if (exist.gdsn(snps_gds_file, "snp_id") == T) 
            read.gdsn(index.gdsn(snps_gds_file, "snp_id")) else 
              read.gdsn(index.gdsn(snps_gds_file, "snp.id"))
          
          rownames(snps_mtrx)   <- if (exist.gdsn(snps_gds_file, "sample_id") == T) 
              read.gdsn(index.gdsn(snps_gds_file, "sample_id")) else 
                read.gdsn(index.gdsn(snps_gds_file, "sample.id"))
          
          closefn.gds(snps_gds_file)
          
          print("Genotype  matrix has been loaded", quote = F)
          print(paste0("Genotype matrix: ", nrow(snps_mtrx), " x ", ncol(snps_mtrx)), quote = F)
          return(snps_mtrx)
        },
        error = function(cond){
          message(cond)
          return(NA)
        }
      )
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
  dnam_mtrx = tryCatch(
    {
      print("Loading DNAm  matrix... ", quote = F)
      
      dnam_gds_file         <- openfn.gds(gds_fn)
      
      dnam_mtrx             <- read.gdsn(index.gdsn(dnam_gds_file, "beta_mtrx"))
      colnames(dnam_mtrx)   <- read.gdsn(index.gdsn(dnam_gds_file, "cpg_id"))
      rownames(dnam_mtrx)   <- read.gdsn(index.gdsn(dnam_gds_file, "sample_id"))
      
      closefn.gds(dnam_gds_file)
      
      print("DNAm  matrix has been loaded", quote = F)
      print(paste0("DNAm matrix: ", nrow(dnam_mtrx), " x ", ncol(dnam_mtrx)), quote = F)
      return(dnam_mtrx)
    },
    error = function(cond){
      print("The DNAm matrix cannot be loaded: ", quote = F)
      message(cond)
      
      return(NA)
    }
  )
  
  return(dnam_mtrx)
}
