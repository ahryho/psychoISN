require(foreign)
require(data.table)

library(parallel)
library(foreach)
library(doParallel)

# Correct for with age, sex, BMI,Smoking Score, DNAm_SV{1-3}, DNA_PC{1,2}

# 1. Load data

# args        <- commandArgs(T)
# mtrx_fn     <- as.character(args[1]) # dnam, snps, gex
# pheno_fn    <- as.character(args[3])
# out_mtrx_fn <- as.character(args[4]) # corrected matrix
# treatment   <- as.character(args[2]) # dex
# 
# mtrx_fn     <- paste0(mtrx_fn, treatment, ".csv")

treatment   <- "veh"

mtrx_fn     <- paste0("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/data/integrative/matrixEQTL/methyl_beta_mtrx_", treatment, ".csv")
pheno_fn    <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/data/pheno/pheno_full_for_kimono.csv"
out_mtrx_fn <- "/binder/mgp/workspace/2020_DexStim_Array_Human/isn/input/dnam/methyl_beta_mtrx_corrected_for_cov"

mtrx    <- fread(mtrx_fn) 
cpg_ids <- mtrx$CpG_ID
mtrx    <- mtrx[, -1]

pheno    <- read.csv2(pheno_fn) #, na.strings = "#N/A") 
pheno    <- pheno[pheno$Include == 1 & pheno$Group == treatment, ]

# 2. Prepare covariates 

cov_lst <- c("Sex", "Age", "BMI_D1", "DNAm_SV1", "DNAm_SV2", "DNAm_SV3", "PC1", "PC2")
cov_df  <- pheno[, c("DNA_ID", cov_lst)]

cov_df  <- cov_df[match(colnames(mtrx), cov_df$DNA_ID ),]
cov_df[, -1]  <- scale(cov_df[, -1]) 
mtrx <- as.matrix(mtrx)

# 3. Build model

no.cores <- detectCores() - 1
cl <- makeCluster(no.cores)
registerDoParallel(cl)

res <- foreach(feature = 1:nrow(mtrx), .combine = rbind) %dopar% { 
  lm.model <- lm(mtrx[feature, ] ~ ., data = cov_df[, -1])
  predict(lm.model)
}

stopImplicitCluster()

res <- as.data.frame(res)
res <- cbind(cpg_ids, mtrx)
colnames(res) <- c("CpG_ID", colnames(mtrx))

# 4. Save results

fwrite(res, 
       paste0(out_mtrx_fn, "_", treatment, ".csv"),
       quote = F, row.names = F, sep = ";")
