library(dplyr)
library(data.table)
library(gdsfmt)

source("~/kul/dex-stim-human-array-isns/code/00_functions/flatten_matrix.R")

dir_pre  <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/"
rslt_dir <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/results/5_fold_cv/all/"
setwd(rslt_dir)

treatment       <- "dex"
dex_global_net  <- readRDS(paste0(rslt_dir, treatment, "/networks/smccnet_global_chr_all.rds"))

treatment       <- "veh"
veh_global_net  <- readRDS(paste0(rslt_dir, treatment, "/networks/smccnet_global_chr_all.rds"))

## Flatten global atrices

dex_flat_df <- flatten_matrix(dex_global_net) %>% dplyr::select(V1 = row, V2 = column, weight_dex = weight)
veh_flat_df <- flatten_matrix(veh_global_net) %>% dplyr::select(V1 = row, V2 = column, weight_veh = weight)

flat_full_df     <- full_join(veh_flat_df, dex_flat_df, by = c("V1", "V2")) %>% 
  replace(is.na(.), 0) %>% setDT()

# flat_inner_df <- flat_full_df[weight_veh & weight_dex != 0] %>% arrange(weight_veh, weight_dex)

flat_df <- rbind(flat_full_df[weight_dex == 0 & weight_veh >= 0.9 ], # disappearing
                 flat_full_df[weight_veh == 0 & weight_dex >= 0.85 ], # appearing
                 flat_full_df[weight_veh >= 0.2 & weight_dex >= 0.2] # remaining
                 )

# fwrite(flat_df, "mnda/global_network_mnda_input.csv", quote = F, sep = ";")

## Prepare ISNs' matrices for MNDA 

associations_of_interest <- flat_df[, 1:2]

### Load a ISNs' filenames: dex and veh have the same names, the difference is in the folder name
isns_fn_lst <- list.files(paste0("dex/networks"), pattern = "individual", full.names = T)

lapply(isns_fn_lst, function(fn){
  sample_id <- gsub(".*individual_(.+).rds", "\\1", fn)
  
  dex_isn  <- readRDS(paste0(rslt_dir, "dex/networks/smccnet_global_chr_all.rds"))
  veh_isn  <- readRDS(paste0(rslt_dir, "veh/networks/smccnet_global_chr_all.rds"))
  
  ## Flatten global atrices
  
  dex_flat_df  <- flatten_matrix(dex_isn) %>% select(V1 = row, V2 = column, weight_dex = weight)
  veh_flat_df  <- flatten_matrix(veh_isn) %>% select(V1 = row, V2 = column, weight_veh = weight)
  
  flat_full_df <- full_join(veh_flat_df, dex_flat_df, by = c("V1", "V2")) %>%
    replace(is.na(.), 0) %>% setDT()
  
  isn_mnda_input <- inner_join(associations_of_interest, flat_full_df)
  
  fwrite(isn_mnda_input, 
         paste0("mnda/isns/isn_mnda_input_", sample_id, ".csv"), 
         quote = F, sep = ";")
})
