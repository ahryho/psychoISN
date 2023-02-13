library(SmCCNet)
library(gdsfmt)

require(data.table)

source("~/kul/dex-stim-human-array-isns/code/00_functions/load_gds_files.R")

setwd("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/")

# 1. Set up global variables

treatment <- "veh"
cv_k      <- 5
chr       <- "all"
cv_dir    <- paste0("results/", cv_k, "_fold_cv/", chr, "/") # chromosomes/", chr)

# dnam_gds_fn <- paste0("input/dnam/mad_filtered/gds/chromosomes/", treatment, "/methyl_beta_mtrx_corrected_for_cov", "_", treatment, "_chr", chr, ".gds")
# snps_gds_fn <- paste0("input/snps/ld_pruned/gds/chromosomes/dex_geno_chr", chr, ".gds")

dnam_gds_fn <- paste0("input/dnam/mad_filtered/gds/methyl_beta_mtrx_corrected_for_cov_mad80_filtered_", treatment, ".gds")
snps_gds_fn <- paste0("input/snps/ld_pruned/gds/dex_geno_imputed_maf_ld_pruned_from_gen.gds")

dnam_mtrx <- LoadMethyl(dnam_gds_fn, is_mad = F)
snps_mtrx <- LoadGenotype(snps_gds_fn, is_ld = F)

pheno           <- fread("input/pheno/pheno_full_for_kimono.csv", dec = ",")[Include == T]
pheno_treatmnet <- pheno[Dex == ifelse(treatment == "veh", 0, 1)]

dnam_mtrx <- dnam_mtrx[pheno_treatmnet$DNA_ID, ]
snps_mtrx <- snps_mtrx[pheno_treatmnet$DNA_ID, ]

# 2. Load SmCCNet results 

modules_obj <- readRDS(file = paste0(cv_dir, treatment, "/smccnet_omic_modules_chr_", chr, ".rds"))
weights     <- modules_obj$weights
sim_mtrx    <- modules_obj$sim_mtrx
modules     <- modules_obj$modules
features    <- modules_obj$non_zero_features
  
nr_cpgs     <- table(features %in% colnames(dnam_mtrx))[2]

# 3. Load genomic coordintaes

## 3.1. Load CpGs' coordinates

cpg_loc <- LoadMethylCoordinates(dnam_gds_fn)

## 3.1. Load SNP' coordinates

snp_loc <- LoadGenotypeCoordinates(snps_gds_fn)

# 4. Check how many cis / trans

features

# 3. Plot modules 

bigCor <- cor(cbind(dnam_mtrx[, features[1:nr_cpgs]], 
                    snps_mtrx[, features[(nr_cpgs + 1):length(features)]]))
edgeCut <- 0.005
for(idx in 1:length(modules)){
  filename <- paste0(cv_dir, chr, "/module_", idx, ".pdf")
  plotMultiOmicsNetwork(Abar = sim_mtrx, CorrMatrix = bigCor,
                        multiOmicsModule = modules, ModuleIdx = idx, P1 = nr_cpgs,
                        EdgeCut = edgeCut, FeatureLabel = features,
                        SaveFile = filename)
}
