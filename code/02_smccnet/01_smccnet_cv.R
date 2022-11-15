setwd("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/")

source("~/kul/dex-stim-human-array-isns/code/00_functions/load_gds_files.R")
source("~/kul/dex-stim-human-array-isns/code/00_functions/load_pkgs.R")

pkgs_list <- c("SmCCNet", "gdsfmt", "SNPRelate",
               "data.table",
               "parallel", "foreach", "doParallel", 
               "tictoc")

LoadPackages(pkgs_list)

# 1. Set up global variables

args        <- commandArgs(T)
treatment   <- as.character(args[1]) 
chr         <- as.numeric(args[2])  # "all_dim_reduction_mad_80"
cv_k        <- as.numeric(args[3])
cv_dir      <- as.character(args[4])
dnam_gds_fn <- as.character(args[5])
snps_gds_fn <- as.character(args[6])

pheno_trait <- c("Status")

# chr <- 22
# treatment <- "veh"
# dnam_gds_fn <- paste0("input/dnam/mad_filtered/gds/chromosomes/", treatment, "/methyl_beta_mtrx_corrected_for_cov", "_", treatment, "_chr", chr, ".gds")
# snps_gds_fn <- paste0("input/snps/ld_pruned/gds/chromosomes/dex_geno_chr", chr, ".gds")
# cv_k    <- 5
# cv_dir  <- paste0("results/", cv_k, "_fold_cv/")
## cv_dir  <- paste0("tmp_data/example_chr_", chr, "_", cv_k, "_fold_cv/")

# 2. Load DNAm beta mtrx

tic("Total computation time of SmCCNet CV")

tic("Load DNAm")
dnam_mtrx <- LoadMethyl(dnam_gds_fn, is_mad = F)
toc()

# 3. Load SNP data

tic("Load Genotype mtrx")
snps_mtrx <- LoadGenotype(snps_gds_fn, is_ld = F)
toc()

# 4. Load pheno trait

pheno           <- fread("input/pheno/pheno_full_for_kimono.csv", dec = ",")[Include == T]
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

cc_coef <- NULL # Unweighted version of SmCCNet.

#### Feature sampling proportions. 
s1      <- 0.9 
s2      <- 0.4 

#### Number of subsamples
subsample_nr <- 50

### 5.2.2. Create sparsity penalty options.
penalty_1 <- seq(.01, .1, by = .01)
penalty_2 <- seq(.01, .1, by = .01)
penalty_grid <- expand.grid(penalty_1, penalty_2)

#### Set a CV directory.
system(paste0("mkdir -p ", cv_dir))

### 5.2.3. Create training and test data

set.seed(1234)

#### For each of the K-fold we compute the prediction err for each penalty pair
res_cv_penaltiy_grid <- matrix(0, nrow = 2, ncol = cv_k * nrow(penalty_grid))

fold_test_idx <- split(1:nr_samples, sample(1:nr_samples, cv_k))

print(paste0("Start ", cv_k, "-fold CV... "), quote = F)
print(paste0("Start date and time: ", Sys.time()), quote = F)
tic("SmCCNet CV")

for(i in 1:cv_k){
  print(paste0("Round ", i, " of CV"))
  
  idx <- fold_test_idx[[i]]
  
  dnam_train <- scale(dnam_mtrx[-idx, ])
  dnam_test  <- scale(dnam_mtrx[idx, ])
  
  snps_train <- scale(snps_mtrx[-idx, ],  scale = apply(snps_mtrx[-idx, ], 2, sd, na.rm = T) + 1e-6)
  snps_test  <- scale(snps_mtrx[idx, ],  scale = apply(snps_mtrx[idx, ], 2, sd, na.rm = T) + 1e-6)
  
  trait_train <- scale(pheno_trait_vec[-idx, ],  scale = apply(pheno_trait_vec[-idx, ], 2, sd, na.rm = T) + 1e-6)
  trait_test <- scale(pheno_trait_vec[idx, ],  scale = apply(pheno_trait_vec[idx, ], 2, sd, na.rm = T) + 1e-6)
  
  # Check if standardized data sets are valid.
  if(is.na(min(min(dnam_train), min(snps_train), min(trait_train), 
               min(dnam_test), min(snps_test), min(trait_test)))){
    stop(paste0("Invalid scaled data. At least one of the data matrices include a
         column with zero variance, iteration = ", idx, "\n"))
  }
  
  no.cores <- as.integer(detectCores() * 0.8)
  cl <- makeCluster(no.cores, type = "FORK")
  registerDoParallel(cl)
  
  clusterEvalQ(cl, library(SmCCNet))
  clusterExport(cl, c("dnam_train", "snps_train", "dnam_test", "snps_test", "trait_train", "trait_test", 
                      "penalty_grid", "s1", "s2", "subsample_nr", "nr_samples", "nr_cpgs", "nr_snps", "cc_coef"))
  
  res <- parSapply(cl, 1:nrow(penalty_grid), function(idx){
    # Consider one pair of sparsity penalties at a time.
    l1 <- penalty_grid[idx, 1]
    l2 <- penalty_grid[idx, 2]
    
    # Run SmCCA on the subsamples (Figure 1, Step II)
    Ws <- getRobustPseudoWeights(dnam_train, snps_train, trait_train, 
                                 l1, l2, s1, s2, 
                                 NoTrait = FALSE,
                                 FilterByTrait = FALSE,
                                 SubsamplingNum = subsample_nr,
                                 CCcoef = cc_coef,
                                 trace = FALSE)
    
    # Aggregate pseudo-canonical weights from the subsamples
    if(is.matrix(Ws)){meanW <- rowMeans(Ws)} else{meanW = Ws}
    v <- meanW[1:nr_cpgs]
    u <- meanW[nr_cpgs + 1:nr_snps]
    
    # Compute the prediction error for given CV fold and sparsity penalties.
    if(is.null(cc_coef)){cc_coef <- rep(1, 3)} # Unweighted SmCCA.

    rho_train <- cor(dnam_train %*% v, snps_train %*% u) * cc_coef[1] +
      cor(dnam_train %*% v, trait_train) * cc_coef[2] +
      cor(snps_train %*% u, trait_train) * cc_coef[3]
    rho_test  <- cor(dnam_test %*% v, snps_test %*% u) * cc_coef[1] +
      cor(dnam_test %*% v, trait_test) * cc_coef[2] +
      cor(snps_test %*% u, trait_test) * cc_coef[3]
    
    delta_cor <- abs(rho_train - rho_test)
    rho_train <- round(rho_train, digits = 5)
    rho_test  <- round(rho_test, digits = 5)
    
    return(list(RhoTest = rho_test, DeltaCor = delta_cor))
  })
  stopCluster(cl)
  
  res_cv_penaltiy_grid[,((i-1)*nrow(penalty_grid)+1):(i*nrow(penalty_grid))] = matrix(unlist(res), nrow = 2)
}

print("Cross-validation has been completed", quote = F)
print(paste0("End date and time: ", Sys.time()), quote = F)
toc()

test_cc    <- matrix(res_cv_penaltiy_grid[1,], nrow = cv_k, byrow = TRUE)
pred_error <- matrix(res_cv_penaltiy_grid[2,], nrow = cv_k, byrow = TRUE)

#### Combine prediction errors from all K folds and compute the total prediction error for each sparsity penalty pair
avg_test_cc     <- colMeans(test_cc)
avg_pred_error  <- colMeans(pred_error)
total_pred_grid <- cbind(penalty_grid, avg_test_cc, avg_pred_error)
colnames(total_pred_grid) <- c("l1", "l2", "test_cc", "pred_error")

print("Saving the total prediction error for each sparsity penalty pair...", quote = F)
tic("Save prediction")

fwrite(total_pred_grid, 
       paste0(cv_dir, "/cv_prediction_grid_", cv_k, "_fold.csv"),
       quote = F, row.names = F, sep = ";")

print(paste("Prediction has been saved into ", cv_dir, "cv_prediction_grid_", cv_k, "_fold.csv"), quote = F)
toc()

#### Extract optimal penalties
min_penalty <- which(avg_pred_error == min(avg_pred_error))
l1 <- total_pred_grid$l1[min_penalty]
l2 <- total_pred_grid$l2[min_penalty]

print(paste0("Optimal penalty pair (lambda 1, lambda 2): (", l1, ", ", l2, ")"))

toc()