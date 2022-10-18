
CreateGDSforSmCCNetCV <- function(gds_fn){
  gds_file <- createfn.gds(gds_fn)
  
  add.gdsn(gds_file, name = "sample_id", val = colnames(mtrx[, -1]))
  
  add.gdsn(gds_file, name = "beta_mtrx", val = t(mtrx[, -1]))
  add.gdsn(gds_file, name = "cpg_id", val = mtrx$CpG_ID)
  
  add.gdsn(gds_file, name = "beta_mtrx", val = t(mtrx[, -1]))
  add.gdsn(gds_file, name = "cpg_id", val = mtrx$CpG_ID)
  
  show(gds_file)
  
  closefn.gds(gds_file)
}