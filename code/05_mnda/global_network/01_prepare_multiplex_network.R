library(dplyr)
library(data.table)
library(gdsfmt)

source("~/kul/dex-stim-human-array-isns/code/00_functions/load_gds_files.R")

dir_pre  <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/"
rslt_dir <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/results/5_fold_cv/all/"
setwd(rslt_dir)

treatment <- "dex"
dex_global_net <- readRDS(paste0(rslt_dir, treatment, "/networks/smccnet_global_chr_all.rds"))

treatment <- "veh"
veh_global_net <- readRDS(paste0(rslt_dir, treatment, "/networks/smccnet_global_chr_all.rds"))
