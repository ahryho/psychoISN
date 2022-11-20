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
chr         <- as.character(args[2])  # "all_dim_reduction_mad_80"
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
# cv_dir  <- paste0("results/", cv_k, "_fold_cv/chromosomes/", chr)
## cv_dir  <- paste0("results/", cv_k, "_fold_cv/")

# 2. Load DNAm beta mtrx

tic("Total computation time of obtaining omic modules")

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

features <- c(colnames(dnam_mtrx), colnames(snps_mtrx))

## 5.2. Set up SmCCNaet parameters

cc_coef  <- NULL # Unweighted version of SmCCNet.

#### Feature sampling proportions. 
s1      <- 0.9 
s2      <- 0.4

#### Number of subsamples
subsample_nr <- 100

#### Extract optimal penalties

total_pred_grid <- fread(paste0(cv_dir, "/cv_prediction_grid_", cv_k, "_fold.csv"))
min_penalty     <- total_pred_grid[pred_error == min(pred_error)]
l1 <- min_penalty$l1[1]
l2 <- min_penalty$l2[1]

print(paste0("Optimal penalty pair (lambda 1, lambda 2): (", l1, ", ", l2, ")"))

# 6. Integrate two omic data types and a quantitative phenotype

print("Start calculation of the canonical correlation weights using optimal penaltiy values ...", quote = F)
print(paste0("Start date and time: ", Sys.time()), quote = F)
tic("Canonical weights computation")

Ws <- getRobustPseudoWeights(dnam_mtrx, snps_mtrx, pheno_trait_vec, l1, l2,
                             s1, s2, NoTrait = FALSE,
                             FilterByTrait = FALSE,
                             SubsamplingNum = subsample_nr,
                             CCcoef = cc_coef)

print("The canonical correlation weights has been calculated", quote = F)
print(paste0("End date and time: ", Sys.time()), quote = F)
toc()

## Extract features which have at least one non-zero value

rownames(Ws) <- features
Ws           <- Ws[apply(Ws, 1, function(x) !all(x==0)),]
new_features <- rownames(Ws)
nr_cpgs      <- table(new_features %in% colnames(dnam_mtrx))[2]

print("Dimenasionality of the similarity matrix after removing zeros: ", quote = F)
print(dim(Ws), quote = F)

print("Saving the canonical correlation (CC) weights...", quote = F)
tic("Save CC weights")

fwrite(Ws,
       paste0(cv_dir, "/smccnet_cc_weights.csv"),
       quote = F, row.names = F, sep = ";")

print(paste("CC weights has been saved into ", cv_dir, "smccnet_cc_weights_chr_", chr, ".csv"), quote = F)
toc()

gc()

# 7. Compute the similarity matrix based on the canonical correlation weight vectors

print("Start computation of the similarity matrix ...", quote = F)
print(paste0("Start date and time: ", Sys.time()), quote = F)
tic("Similarity matrix computation")

sim_mtrx <- getAbar(Ws, FeatureLabel = new_features)

print("The similarity matrix has been calculated", quote = F)
print(paste0("End date and time: ", Sys.time()), quote = F)
toc()

gc()

# 8.  Obtain multi-omics modules

print("Start obtaining multi-omics modules ...", quote = F)
print(paste0("Start date and time: ", Sys.time()), quote = F)
tic("Obtain mult-mics modules")

modules <- getMultiOmicsModules(sim_mtrx, nr_cpgs)

print("Multi-omics modules have been obtained", quote = F)
print(paste0("End date and time: ", Sys.time()), quote = F)
toc()

print("Start saving multi-omics modules ...", quote = F)
tic("Save mult-mics modules as RDS object")

saveRDS(list(weights = Ws, non_zero_features = new_features, sim_mtrx = sim_mtrx, modules = modules), 
        file = paste0(cv_dir, "/smccnet_omic_modules_chr_", chr, ".rds"))

print(paste0("Multi-omics modules have been saved into ", cv_dir), quote = F)
print(paste0("End date and time: ", Sys.time()), quote = F)
toc()

toc()