library(SmCCNet)
library(gdsfmt)
library(SNPRelate)

setwd("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/")

source("~/kul/dex-stim-human-array-isns/code/00_functions/load_gds_files.R")

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
snps_mtrx    <- LoadGenotype(snps_gds_fn, is_ld = F)

# 4. Load pheno trait

pheno           <- fread("input/pheno/pheno_full_for_kimono.csv")
pheno           <- pheno[Include == T]
pheno_treatmnet <- pheno[Dex == ifelse(treatment == "veh", 0, 1)]
pheno_trait_vec <- pheno_treatmnet[, ..pheno_trait]

# 5. SmCCNet

## 5.1. Set up variables

nr_cpgs    <- ncol(dnam_mtrx)
nr_snps    <- ncol(snps_mtrx)
nr_samples <- nrow(snps_mtrx)

dnam_mtrx <- dnam_mtrx[pheno_treatmnet$DNA_ID, ]
snps_mtrx <- snps_mtrx[pheno_treatmnet$DNA_ID, ]

features  <- c(colnames(dnam_mtrx), colnames(snps_mtrx))

## 5.2 Determine optimal sparsity penalties through CV 

### 5.2.1. Set up parameters

cv_k <- 5
cc_coef <- NULL # Unweighted version of SmCCNet.
s1 <- 0.7; s2 <- 0.9 # Feature sampling proportions. 
subsample_nr <- 100 # Number of subsamples.

# Create sparsity penalty options.
pen1 <- seq(.05, .3, by = .05); pen2 <- seq(.05, .3, by = .05) 
P1P2 <- expand.grid(pen1, pen2) # Map (l1, l2) to (c1, c2).
c1 <- sqrt(p1 * s1) * P1P2[ , 1]; c1[c1] <- 1
c2 <- sqrt(p2 * s2) * P1P2[ , 2]; c2[c2 < 1] <- 1
# Based on prior knowledge we may assume that there are at least as many genes # as miRNAs in each network.
P1P2 <- P1P2[which(c1>c2), ]
# Set a CV directory.
cv_dir <- paste0("tmp_data/example_", cv_k, "_fold_cv/")
dir.create(cv_dir)

### 5.2.2. Create training and test data

set.seed(4828)

fold_test_idx <- split(1:nr_samples, sample(1:nr_samples, cv_k))

for(i in 1:cv_k){
  idx <- fold_idx[[i]]
  
  dnam_train <- scale(dnam_mtrx[-idx, ])
  dnam_test  <- scale(dnam_mtrx[idx, ])
  
  snps_train <- scale(snps_mtrx[-idx, ])
  snps_test  <- scale(snps_mtrx[idx, ])
  
  trait_train <- scale(pheno_trait_vec[-idx, ])
  trait_test <- scale(pheno_trait_vec[idx, ])
  
  # Check if standardized data sets are valid.
  if(is.na(min(min(dnam_train), min(snps_train), min(trait_train), 
               min(dnam_test), min(snps_test), min(trait_test)))){
    stop("Invalid scaled data. At least one of the data matrices include a
         column with zero variance.")
  }
  
  sub_dir <- paste0(cv_dir, "cv_", i, "th")
  dir.create(sub_dir)
  
  save(x1.train, x2.train, yy.train, x1.test, x2.test, yy.test,
       s1, s2, P1P2, p1, p2, SubsamplingNum, CCcoef,
       file = paste0(subD, "Data.RData"))
}


l1 <- P1P2[idx, 1]
l2 <- P1P2[idx, 2]
Ws <- getRobustPseudoWeights(dnam_train[, 1:10], snps_train[,1:30], trait_train, l1, l2,
                             s1, s2, NoTrait = FALSE,
                             FilterByTrait = FALSE,
                             SubsamplingNum = subsample_nr,
                             CCcoef = cc_coef)
