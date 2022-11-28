library(dplyr)
library(data.table)
library(gdsfmt)

source("~/kul/dex-stim-human-array-isns/code/00_functions/load_gds_files.R")
source("~/kul/dex-stim-human-array-isns/code/00_functions/get_chord_diagramm.R")

rslt_dir <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/results/5_fold_cv/all/"
setwd(rslt_dir)

# 1. Load global network

treatment <- "dex"
dex_global_net <- readRDS(paste0(treatment, "/networks/smccnet_global_chr_all.rds"))

treatment <- "veh"
veh_global_net <- readRDS(paste0(treatment, "/networks/smccnet_global_chr_all.rds"))

# 2. Check the overlap between dex and veh features

dex_features <- rownames(dex_global_net) # 5'388
veh_features <- rownames(veh_global_net) # 4'519

olap <- intersect(dex_features, veh_features) # 758
olap_cpgs <- olap[grepl("cg", olap, ignore.case = T)] # 566
olap_snps <- olap[grepl("rs", olap, ignore.case = T)] # 192

# 3. Load genomic coordintaes

## 3.1. Load CpGs' coordinates
setwd("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/")

treatment   <- "veh"
dnam_gds_fn <- paste0("input/dnam/mad_filtered/gds/methyl_beta_mtrx_corrected_for_cov_mad80_filtered_", treatment, ".gds")
cpg_loc     <- LoadMethylCoordinates(dnam_gds_fn) %>% data.frame() %>% setDT()

## 3.1. Load SNP' coordinates

snps_gds_fn <- paste0("input/snps/ld_pruned/gds/dex_geno_imputed_maf_ld_pruned_from_gen.gds")
snp_loc     <- LoadGenotypeCoordinates(snps_gds_fn) %>% data.frame() %>% setDT()


####------------
# 4. Check how many cis / trans

# modules_features <- features[do.call(c, modules)]
# 
# cpg_loc[CpG_ID %in% modules_features]
# snp_loc[SNP %in% modules_features]

####------------

# 5. Coord plot 

## 5.1. Baseline

net <- veh_global_net %>% as.matrix() %>% as.data.frame()
net <- veh_global_net[!rownames(veh_global_net) %in% olap,
                      !colnames(veh_global_net) %in% olap] %>% as.matrix() %>% as.data.frame()
CustomChordDiagram(net, cpg_loc, snp_loc, thrsh = 0.80)

## 5.2. Post-dexamethasone

net <- dex_global_net %>% as.matrix() %>% as.data.frame()
CustomChordDiagram(net, cpg_loc, snp_loc, thrsh = 0.25)

## 5.3. Net with olap features betweeb pre- and post-dex

### 5.3.1 Baseline

net <- veh_global_net[olap, olap] %>% as.matrix() %>% as.data.frame()
CustomChordDiagram(net, cpg_loc, snp_loc, thrsh = 0)





# Get CpGs' and SNPs' locations for only features present in the network "net"
net_cpg_loc <- cpg_loc[CpG_ID %in% colnames(net)] %>% mutate(ID = CpG_ID) %>% select(-CpG_ID)
net_snp_loc <- snp_loc[SNP %in% colnames(net)] %>% mutate(ID = SNP) %>% select(-SNP)
net_loc     <- rbind(net_cpg_loc, net_snp_loc)[match(rownames(net), ID)]
