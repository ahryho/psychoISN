# 1. Transform DNAm and SNPs data: individuals are in the rows, features are in the columns

## Load libraries

library(data.table)

## Load data

dnam_base_mtrx <- fread("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/data/integrative/matrixEQTL/methyl_beta_mtrx_veh.csv")
dnam_dex_mtrx  <- fread("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/data/integrative/matrixEQTL/methyl_beta_mtrx_dex.csv")

## Transform the data: individuals are in the rows, features are in the columns

dnam_base_mtrx_t <- t(dnam_base_mtrx)
dnam_dex_mtrx_t  <- t(dnam_dex_mtrx)

## Save matrices
fwrite(dnam_base_mtrx_t,
       "/binder/mgp/datasets/2020_DexStim_Array_Human/isn/methyl_beta_mtrx_veh.csv",
       quote = F, row.names = F, sep = ";")

fwrite(dnam_dex_mtrx_t,
       quote = F, row.names = F, sep = ";")

### The data will be stored in /binder/mgp/datasets/2020_DexStim_Array_Human/isn

# Load the data from dex-stim-human-array proejct

# 2. Correct for covariates