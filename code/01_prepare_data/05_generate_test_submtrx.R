library(SmCCNet)
library(gdsfmt)
library(SNPRelate)

setwd("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns")

source("~/kul/dex-stim-human-array-isns/code/util.R")

# 1. Set up global variables

treatment   <- "veh"

# 2. Load DNAm beta mtrx

gds_fn    <- paste0("input/dnam/gds/methyl_beta_mtrx_corrected_for_cov",  
                    "_", treatment, ".gds") 

dnam_mtrx <- LoadMethyl(gds_fn, is_mad = F)

# 3. Generate test DNAm mtrx

nr_all_cpgs     <- ncol(dnam_mtrx)
nr_desired_cpgs <- 100

dnam_sub_mtrx <- dnam_mtrx[, sample(x = 1:nr_all_cpgs, size = nr_desired_cpgs, replace = F)]

## Generate GDS object and save

gds_fn   <- paste0("input/test_data/methyl_beta_mtrx_corrected_for_cov",  
                   "_", treatment, ".gds") 

gds_file <- createfn.gds(gds_fn)

add.gdsn(gds_file, name = "beta_mtrx", val = dnam_sub_mtrx)
add.gdsn(gds_file, name = "cpg_id", val = colnames(dnam_sub_mtrx))
add.gdsn(gds_file, name = "sample_id", val = rownames(dnam_sub_mtrx))

show(gds_file)

closefn.gds(gds_file)

# 4. Load SNP data

snps_gds_fn <- "input/snps/gds/dex_geno_imputed.gds"
snps_mtrx    <- LoadGenotype(snps_gds_fn, is_ld = F)

# 5. Generate test SNP mtrx

nr_all_snps     <- ncol(snps_mtrx)
nr_desired_snps <- 300

snps_sub_mtrx <- snps_mtrx[, sample(x = 1:nr_all_snps, size = nr_desired_snps, replace = F)]

## Generate GDS object and save

gds_fn   <- "input/test_data/dex_geno_imputed.gds"

gds_file <- createfn.gds(gds_fn)

add.gdsn(gds_file, name = "genotype", val = snps_sub_mtrx)
add.gdsn(gds_file, name = "snp_id", val = colnames(snps_sub_mtrx))
add.gdsn(gds_file, name = "sample_id", val = rownames(snps_sub_mtrx))

show(gds_file)

closefn.gds(gds_file)
