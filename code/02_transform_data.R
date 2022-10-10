# 1. Transform DNAm and SNPs data: individuals are in the rows, features are in the columns
# The data will be stored in /binder/mgp/workspace/2020_DexStim_Array_Human/isn

## Load libraries

library(data.table)
library(arrow)
library(dplyr)

## Function to transform data

transform_data <- function(mtrx, col_id, out_fn = ""){
  # column "col_id" of mtrx is ID
  mtrx_t <- data.frame(t(select(mtrx, -all_of(col_id))))
  colnames(mtrx_t) <- as.data.frame(mtrx)[, col_id]
  mtrx_t[["DNA_ID"]] <- rownames(mtrx_t)
  
  # write_csv_arrow(mtrx_t, out_fn, 
  #                  include_header = T)
  
  # write_parquet(mtrx_t, out_fn)
  
  return(mtrx_t)
}

## Baseline
dnam_base_mtrx <- read_delim_arrow("/binder/mgp/workspace/2020_DexStim_Array_Human/isn/input/dnam/methyl_beta_mtrx_corrected_for_cov_veh.csv", delim = ";")
                                   # /binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/data/integrative/matrixEQTL/methyl_beta_mtrx_veh.csv")

out_fn <- "/binder/mgp/workspace/2020_DexStim_Array_Human/isn/input/dnam/methyl_beta_mtrx_veh"

mtrx_base <- transform_data(dnam_base_mtrx, "CpG_ID", paste0(out.fn, ".parquet")) 

write_feather(mtrx_base, paste0(out_fn, ".arrow"))

## Dex
dnam_dex_mtrx  <-  read_delim_arrow("/binder/mgp/workspace/2020_DexStim_Array_Human/isn/input/dnam/methyl_beta_mtrx_corrected_for_cov_dex.csv", delim = ";")

mtrx_dex <- transform_data(dnam_dex_mtrx, "CpG_ID", 
                       "/binder/mgp/workspace/2020_DexStim_Array_Human/isn/input/dnam/methyl_beta_mtrx_dex.csv") 
