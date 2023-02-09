library(dplyr)
library(data.table)
library(rdist)

source("~/kul/dex-stim-human-array-isns/code/00_functions/load_gds_files.R")
source("~/kul/dex-stim-human-array-isns/code/00_functions/save_isns_as_single_object.R")
source("~/kul/dex-stim-human-array-isns/code/00_functions/euclidean_dist.R")

# 1. Set up global variables

args        <- commandArgs(T)
treatment   <- as.character(args[1]) 
rslt_dir    <- as.character(args[2]) 

# rslt_dir <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/results/5_fold_cv/all/"
# treatment <- "veh"
setwd(rslt_dir)

# # 1. Load global network
# 
# treatment <- "dex"
# dex_global_net  <- readRDS(paste0(treatment, "/networks/smccnet_global_chr_all.rds"))
# 
# treatment <- "veh"
# veh_global_net <- readRDS(paste0(treatment,  "/networks/smccnet_global_chr_all.rds"))
# 
# dex_global_net_flatten <- flattenCorrMatrix(dex_global_net, pmat = NULL) %>% setDT()
# veh_global_net_flatten <- flattenCorrMatrix(veh_global_net, pmat = NULL) %>% setDT()
# 
# # 2. Get the overlap between dex and veh features
# 
# dex_features <- rownames(dex_global_net) # 5'388
# veh_features <- rownames(veh_global_net) # 4'519
# 
# olap_feat <- intersect(dex_features, veh_features) # 758
# olap_cpgs <- olap_feat[grepl("cg", olap_feat, ignore.case = T)] # 566
# olap_snps <- olap_feat[grepl("rs", olap_feat, ignore.case = T)] # 192

# 2. Load ISNs

isns_fn_lst <- list.files(paste0(treatment, "/networks"), pattern = "individual", full.names = T)
isns_lst    <- lapply(isns_fn_lst, readRDS)

isns_membership_df <- lapply(isns_fn_lst, function(fn) {
  sample_id <- gsub(".*individual_(.+)_chr_all.rds", "\\1", fn) })

isns_graph_lst     <- isns_lst

# isns_membership_df <- lapply(isns_lst, function(net_obj) data.frame(netID = net_obj[[1]], treatment = treatment)) %>% 
#   bind_rows()
# isns_graph_lst     <- lapply(isns_lst, function(net_obj) as.matrix(net_obj[[2]]))

rm(isns_lst)
gc()

# 3. Calculate distances

isns_dist_mtrx  <- euclidean_dist(G = isns_graph_lst, meth = "euclidean")
colnames(isns_dist_mtrx) <- isns_membership_df #$netID

# 4. Save results

fwrite(as.data.table(isns_dist_mtrx, keep.rownames = F), 
       paste0(rslt_dir, treatment, "/smccnet_isns_euclidean_dist_mtrx.csv"),
       quote = F, sep = ";", row.names = F, col.names = T)
