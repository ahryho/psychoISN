library(arrow)
library(gdsfmt)

# 1. Load corrected for covariates dnam matrices

treatment   <- "dex"

mtrx_fn_prefix <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/input/dnam/methyl_beta_mtrx_corrected_for_cov"
mtrx_fn        <-  paste0(mtrx_fn_prefix, "_", treatment, ".csv")

mtrx <- read_delim_arrow(mtrx_fn,  delim = ";")

# 2. Create a GDS object

gds_fn   <- paste0("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/input/dnam/gds/methyl_beta_mtrx_corrected_for_cov",  
                   "_", treatment, ".gds") 

gds_file <- createfn.gds(gds_fn)

add.gdsn(gds_file, name = "beta_mtrx", val = t(mtrx[, -1]))
add.gdsn(gds_file, name = "cpg_id", val = mtrx$CpG_ID)
add.gdsn(gds_file, name = "sample_id", val = colnames(mtrx[, -1]))

show(gds_file)

closefn.gds(gds_file)
