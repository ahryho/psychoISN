library(SmCCNet)
library(gdsfmt)
library(SNPRelate)

setwd("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/")

source("~/kul/dex-stim-human-array-isns/code/util.R")

# 1. Set up global variables

treatment   <- "veh"
pheno_trait <- "Status"

# 2. Load DNAm beta mtrx
# 
# gds_fn    <- paste0("input/dnam/gds/methyl_beta_mtrx_corrected_for_cov",  
#                     "_", treatment, ".gds") 

gds_fn    <- paste0("input/test_data/methyl_beta_mtrx_corrected_for_cov",  
                    "_", treatment, ".gds") 

dnam_mtrx <- LoadMethyl(gds_fn, is_mad = F)

# 3. Load SNP data

# snps_gds_fn <- "input/snps/gds/dex_geno_imputed.gds"
snps_gds_fn  <- "input/test_data/dex_geno_imputed.gds"
snp_mtrx     <- LoadGenotype(snps_gds_fn, is_ld = F)

# 4. Load pheno trait

pheno           <- fread("input/pheno/pheno_full_for_kimono.csv")
pheno           <- pheno[Include == T]
pheno_treatmnet <- pheno[Dex == ifelse(treatment == "veh", 0, 1)]
pheno_trait_vec <- pheno_treatmnet[, ..pheno_trait]

# 5. SmCCNet

