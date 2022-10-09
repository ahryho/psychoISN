# 1. Transform DNAm and SNPs data: individuals are in the rows, features are in the columns
# The data will be stored in /binder/mgp/datasets/2020_DexStim_Array_Human/isn

## Load libraries

library(data.table)

## Function to transform data

transform_data <- function(mtrx, col_id, out_fn){
  # column "col_id" of mtrx is ID
  mtrx_t <- data.frame(t(select(mtrx, -col_id)))
  colnames(mtrx_t) <- as.data.frame(mtrx)[, col_id]
  
  fwrite(mtrx_t, out_fn, 
         quote = F, row.names = T, sep = ";")
  
  return(mtrx_t)
}

## Baseline
dnam_base_mtrx <- fread("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/data/integrative/matrixEQTL/methyl_beta_mtrx_veh.csv")

mtrx <- transform_data(dnam_base_mtrx, "CpG_ID", 
                       "/binder/mgp/datasets/2020_DexStim_Array_Human/isn/input/dnam/methyl_beta_mtrx_veh.csv") 

## Dex
dnam_dex_mtrx  <- fread("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/data/integrative/matrixEQTL/methyl_beta_mtrx_dex.csv")

mtrx_dex <- transform_data(dnam_dex_mtrx, "CpG_ID", 
                       "/binder/mgp/datasets/2020_DexStim_Array_Human/isn/input/dnam/methyl_beta_mtrx_dex.csv") 
