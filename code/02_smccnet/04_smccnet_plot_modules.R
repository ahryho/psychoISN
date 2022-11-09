library(SmCCNet)
library(gdsfmt)
library(SNPRelate)

require(data.table)

library(parallel)
library(foreach)
library(doParallel)

source("~/kul/dex-stim-human-array-isns/code/00_functions/load_gds_files.R")

setwd("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/")

# 1. Set up global variables

treatment <- "veh"
cv_k      <- 5
cv_dir    <- paste0("results/", cv_k, "_fold_cv/chromosomes/")
chr       <- 1

dnam_gds_fn <- paste0("input/dnam/mad_filtered/gds/chromosomes/", treatment, "/methyl_beta_mtrx_corrected_for_cov", "_", treatment, "_chr", chr, ".gds")
snps_gds_fn <- paste0("input/snps/ld_pruned/gds/chromosomes/dex_geno_chr", chr, ".gds")

dnam_mtrx <- LoadMethyl(dnam_gds_fn, is_mad = F)
snps_mtrx <- LoadGenotype(snps_gds_fn, is_ld = F)

pheno           <- fread("input/pheno/pheno_full_for_kimono.csv", dec = ",")[Include == T]
pheno_treatmnet <- pheno[Dex == ifelse(treatment == "veh", 0, 1)]

dnam_mtrx <- dnam_mtrx[pheno_treatmnet$DNA_ID, ]
snps_mtrx <- snps_mtrx[pheno_treatmnet$DNA_ID, ]

nr_cpgs   <- ncol(dnam_mtrx)
features  <- c(colnames(dnam_mtrx), colnames(snps_mtrx))

# 2. Load SmCCNet results 

modules_obj <- readRDS(file = paste0(cv_dir, chr, "/smccnet_omic_modules_chr_", chr, ".rds"))
sim_mtrx    <- modules_obj$sim_mtrx
modules     <- modules_obj$modules

# 3. Plot modules 

bigCor <- cor(cbind(dnam_mtrx, snps_mtrx))
edgeCut <- 0.005
for(idx in 1:length(modules)){
  filename <- paste0(cv_dir, chr, "/module_", idx, ".pdf")
  plotMultiOmicsNetwork(Abar = sim_mtrx, CorrMatrix = bigCor,
                        multiOmicsModule = modules, ModuleIdx = idx, P1 = nr_cpgs,
                        EdgeCut = edgeCut, FeatureLabel = features,
                        SaveFile = filename)
}
