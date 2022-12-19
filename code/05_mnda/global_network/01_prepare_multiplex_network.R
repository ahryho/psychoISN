library(dplyr)
library(data.table)
library(gdsfmt)

source("~/kul/dex-stim-human-array-isns/code/00_functions/flatten_matrix.R")

dir_pre  <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/"
rslt_dir <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/results/5_fold_cv/all/"
setwd(rslt_dir)

treatment <- "dex"
dex_global_net <- readRDS(paste0(rslt_dir, treatment, "/networks/smccnet_global_chr_all.rds"))

treatment <- "veh"
veh_global_net <- readRDS(paste0(rslt_dir, treatment, "/networks/smccnet_global_chr_all.rds"))

## Flatten matrices

dex_flat_df <- flatten_matrix(dex_global_net)
veh_flat_df <- flatten_matrix(veh_global_net)

flat_full_df     <- full_join(veh_flat_df, dex_flat_df, by = c("row", "column")) %>% 
  select(V1 = row, V2 = column, weight_veh = weight.x, weight_dex = weight.y) %>%
  replace(is.na(.), 0) %>% setDT()

# flat_inner_df <- flat_full_df[weight_veh & weight_dex != 0] %>% arrange(weight_veh, weight_dex)


flat_df <- rbind(flat_full_df[weight_dex == 0 & weight_veh >= 0.85 ], # disappearing
                 flat_full_df[weight_veh == 0 & weight_dex >= 0.85 ], # appearing
                 flat_full_df[weight_veh >= 0.1 & weight_dex >= 0.1] # remaining
                 )

fwrite(flat_df, "mnda/global_network/global_network_mnda_input_2K.csv", quote = F, sep = ";")
