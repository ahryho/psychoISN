library(dplyr)
library(data.table)
library(gdsfmt)
library(RColorBrewer)

library(annotatr)

source("~/kul/dex-stim-human-array-isns/code/00_functions/load_gds_files.R")
source("~/kul/dex-stim-human-array-isns/code/00_functions/annotate_ntwrk_features.R")

dir_pre  <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/"
rslt_dir <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/results/5_fold_cv/all/"
setwd(rslt_dir)

## Load annotation files

cpg_chipseeker_anno <- readRDS(paste0(dir_pre, "input/annotation/cpgs/cpg_annotated_withChIPseeker_all.rds"))
cpg_chipseeker_anno <- cpg_chipseeker_anno@anno

snp_chipseeker_anno <- readRDS(paste0(dir_pre, "input/annotation/snps/snp_annotated_withChIPseeker_all.rds"))
snp_chipseeker_anno <- snp_chipseeker_anno@anno

chromhmm_blood_anno <- readRDS(paste0(dir_pre, "input/annotation/chromHMM/chromHMM_blood_states.Rds"))

## Annotate dex network features

treatment <- "dex"
net            <- readRDS(paste0(rslt_dir, treatment, "/networks/smccnet_global_chr_all.rds"))

### CpGs

out_fn         <- paste0(rslt_dir, treatment, "/smccnet_", treatment, "_cpgs_annotated.rds")
dex_cpg_anno   <- annotate_ntwk_cpgs(net, cpg_chipseeker_anno, out_fn)

## ChromHMM annotation

out_fn                <- paste0(rslt_dir, treatment, "/smccnet_", treatment, "_cpgs_chromhmm_annotated.rds")
dex_cpg_chromhmm_anno <- annotate_chromHMM(dex_cpg_anno, chromhmm_blood_anno, out_fn)

### SNPs

out_fn         <- paste0(rslt_dir, treatment, "/smccnet_", treatment, "_snps_annotated.rds")
dex_snp_anno   <- annotate_ntwk_snps(net, snp_chipseeker_anno, out_fn)

fwrite(list(unique(dex_snp_anno$SYMBOL)), 
       paste0(rslt_dir, treatment, "/smccnet_", treatment, "_snp_gene_symbol_list.csv"),
       sep = ";", quote = F, row.names = F, col.names = F)

## ChromHMM annotation

out_fn                <- paste0(rslt_dir, treatment, "/smccnet_", treatment, "_snps_chromhmm_annotated.rds")
dex_snp_chromhmm_anno <- annotate_chromHMM(dex_snp_anno, chromhmm_blood_anno, out_fn)

## Annotate baseline network features

treatment <- "veh"
net            <- readRDS(paste0(rslt_dir, treatment, "/networks/smccnet_global_chr_all.rds"))

### CpGs

out_fn         <- paste0(rslt_dir, treatment, "/smccnet_", treatment, "_cpgs_annotated.rds")
veh_cpg_anno   <- annotate_ntwk_cpgs(net, cpg_chipseeker_anno, out_fn)

## ChromHMM annotation

out_fn                <- paste0(rslt_dir, treatment, "/smccnet_", treatment, "_cpgs_chromhmm_annotated.rds")
veh_cpg_chromhmm_anno <- annotate_chromHMM(veh_cpg_anno, chromhmm_blood_anno, out_fn)

### SNPs

out_fn         <- paste0(rslt_dir, treatment, "/smccnet_", treatment, "_snps_annotated.rds")
veh_snp_anno   <- annotate_ntwk_snps(net, snp_chipseeker_anno, out_fn)

## ChromHMM annotation

out_fn                <- paste0(rslt_dir, treatment, "/smccnet_", treatment, "_snps_chromhmm_annotated.rds")
veh_snp_chromhmm_anno <- annotate_chromHMM(veh_snp_anno, chromhmm_blood_anno, out_fn)


## Annotate ad save bkgr snps

snps_gds_fn <- paste0("input/snps/ld_pruned/gds/dex_geno_imputed_maf_ld_pruned_from_gen.gds")
snps_mtrx   <- LoadGenotype(snps_gds_fn, is_ld = F)
bkgr_snps   <- snp_chipseeker_anno[names(snp_chipseeker_anno) %in% colnames(snps_mtrx),]

fwrite(list(unique(bkgr_snps$SYMBOL)), 
       paste0(rslt_dir, "/smccnet_bkgr_snp_gene_symbol_list.csv"),
       sep = ";", quote = F, row.names = F, col.names = F)
