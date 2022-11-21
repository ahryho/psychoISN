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

cpg_loc <- LoadMethylCoordinates(dnam_gds_fn) %>% data.frame() %>% setDT()

## 3.1. Load SNP' coordinates

snp_loc <- LoadGenotypeCoordinates(snps_gds_fn) %>% data.frame() %>% setDT()

# 4. Check how many cis / trans

modules_features <- features[do.call(c, modules)]

modules_features <- features[modules[[7]]]

cpg_loc[CpG_ID %in% modules_features]
snp_loc[SNP %in% modules_features]
