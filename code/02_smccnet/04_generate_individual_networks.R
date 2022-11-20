source("~/kul/dex-stim-human-array-isns/code/00_functions/load_pkgs.R")
source("~/kul/dex-stim-human-array-isns/code/00_functions/load_gds_files.R")

pkgs_list <- c("SmCCNet", "gdsfmt", 
               "data.table", "dplyr",
               "parallel", "foreach", "doParallel", 
               "tictoc")

LoadPackages(pkgs_list)

setwd("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/")

# 1. Set up global variables

# 1. Set up global variables

args        <- commandArgs(T)
treatment   <- as.character(args[1]) 
chr         <- as.character(args[2]) 
cv_dir      <- as.character(args[3])
dnam_gds_fn <- as.character(args[4])
snps_gds_fn <- as.character(args[5])

pheno_trait <- c("Status")
# treatment   <- "veh"
# chr         <- "all"
# cv_dir      <- paste0("results/5_fold_cv/all/veh/")
# dnam_gds_fn <- paste0("input/dnam/mad_filtered/gds/methyl_beta_mtrx_corrected_for_cov_mad80_filtered_veh.gds")
# snps_gds_fn <- paste0("input/snps/ld_pruned/gds/dex_geno_imputed_maf_ld_pruned_from_gen.gds")

dnam_mtrx   <- LoadMethyl(dnam_gds_fn, is_mad = F)
snps_mtrx   <- LoadGenotype(snps_gds_fn, is_ld = F)

pheno           <- fread("input/pheno/pheno_full_for_kimono.csv", dec = ",")[Include == T]
pheno_treatmnet <- pheno[Dex == ifelse(treatment == "veh", 0, 1)]
pheno_trait_vec <- pheno_treatmnet[, ..pheno_trait]

dnam_mtrx       <- dnam_mtrx[pheno_treatmnet$DNA_ID, ]
snps_mtrx       <- snps_mtrx[pheno_treatmnet$DNA_ID, ]

# 2. Load SmCCNet results 

modules_obj <- readRDS(file = paste0(cv_dir, "/smccnet_omic_modules_chr_", chr, ".rds"))
weights     <- modules_obj$weights
sim_mtrx    <- modules_obj$sim_mtrx
modules     <- modules_obj$modules
features    <- modules_obj$non_zero_features

# 3. Extract the features obtained by SmCCNet

modules_features <- features[do.call(c, modules)]

# 4. Subset DNAm and SNPs matrices

global_dnam_mtrx <- dnam_mtrx[, colnames(dnam_mtrx) %in% modules_features] 
global_snps_mtrx <- snps_mtrx[, colnames(snps_mtrx) %in% modules_features] 

# 5. Set up SmCCNet parameters

## Sparsity parameter: l1 = l2 : no sparsity, just CCA

l1 <- l2 <- 1

## Sampling proportions: s1 <- s2 <- 1 

s1 <- s2 <- 1

## Subsampling nr: 1

subsample_nr <- 1

# 6. Calculatethe canonical correlation weights

print("Start calculation of the canonical correlation weights ...", quote = F)
print(paste0("Start date and time: ", Sys.time()), quote = F)
tic("Canonical weights computation")

Ws <- getRobustPseudoWeights(global_dnam_mtrx, global_snps_mtrx, pheno_trait_vec, l1, l2,
                             s1, s2, NoTrait = FALSE,
                             FilterByTrait = FALSE,
                             SubsamplingNum = subsample_nr, trace = T)

print("The canonical correlation weights has been calculated", quote = F)
print(paste0("End date and time: ", Sys.time()), quote = F)
toc()

# 7. Compute the similarity matrix based on the canonical correlation weight vectors

print("Start computation of the global similarity matrix ...", quote = F)
print(paste0("Start date and time: ", Sys.time()), quote = F)
tic("Global similarity matrix computation")

global_sim_mtrx <- getAbar(Ws, FeatureLabel = modules_features)

print("The global similarity matrix has been calculated", quote = F)
print(paste0("End date and time: ", Sys.time()), quote = F)
toc()

system(paste0("mkdir -p ", cv_dir, "/networks"))

saveRDS(global_sim_mtrx, 
        file = paste0(cv_dir, "/networks/smccnet_global_chr_", chr, ".rds"))

# 8. Generate individual netwroks

no.cores <- as.integer(detectCores() * 0.8)
cl <- makeCluster(no.cores, type = "FORK")
registerDoParallel(cl)

clusterEvalQ(cl, library(SmCCNet))
clusterExport(cl, c("global_dnam_mtrx", "global_snps_mtrx", "pheno_trait_vec",  
                    "l1", "l2", "s1", "s2", "subsample_nr", "modules_features"))

networks <- parSapply(cl, 1:nrow(pheno_trait_vec), function(sample_idx){
  sample_dna_id <- rownames(global_dnam_mtrx)[sample_idx]
  
  print(paste0("Obtaining an individual network for ", sample_dna_id, " ..."), quote = F)
  
  minus_sample_dnam_mtrx      <- global_dnam_mtrx[-sample_idx, ]
  minus_sample_snps_mtrx      <- global_snps_mtrx[-sample_idx, ]
  minus_samplepheno_trait_vec <- pheno_trait_vec[-sample_idx,]
  
  print("Calculation of the canonical correlation weights ...", quote = F)
  Ws <- getRobustPseudoWeights(minus_sample_dnam_mtrx, minus_sample_snps_mtrx, minus_samplepheno_trait_vec, l1, l2,
                               s1, s2, NoTrait = FALSE,
                               FilterByTrait = FALSE,
                               SubsamplingNum = subsample_nr, trace = T)
  
  print("Computation of the similarity matrix ...", quote = F)
  minus_sample_sim_mtrx <- getAbar(Ws, FeatureLabel = modules_features)
  
  print("Creating an individual network ", quote = F)
  sample_sim_mtrx <- global_sim_mtrx - minus_sample_sim_mtrx
  
  print("Saving an individual network ", quote = F)
  saveRDS(list(DNA_ID = sample_dna_id, network = sample_sim_mtrx), 
          file = paste0(cv_dir, "/networks/smccnet_individual_", sample_dna_id, "_chr_", chr, ".rds"))
  
  return(sample_sim_mtrx)
})

stopCluster(cl)
