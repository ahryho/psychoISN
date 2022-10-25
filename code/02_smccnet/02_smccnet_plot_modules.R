library(SmCCNet)
library(gdsfmt)
library(SNPRelate)

require(data.table)

library(parallel)
library(foreach)
library(doParallel)

setwd("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/")

# 1. Set up global variables

treatment <- "veh"
cv_k      <- 5
cv_dir    <- paste0("tmp_data/example_", cv_k, "_fold_cv/")

# 2. Load SmCCNet results 

modules <- readRDS(file = paste0(cv_dir, "smccnet_res_", cv_k, ".rds"))

# 3. Plot modules 

bigCor <- cor(cbind(dnam_mtrx, snps_mtrx))
edgeCut <- 0.005
for(idx in 1:length(Modules)){
  filename <- paste0(cv_dir, "module_", idx, ".pdf")
  plotMultiOmicsNetwork(Abar = sim_mtrx, CorrMatrix = bigCor,
                        multiOmicsModule = Modules, ModuleIdx = idx, P1 = nr_cpgs,
                        EdgeCut = edgeCut, FeatureLabel = features,
                        SaveFile = filename)
}
