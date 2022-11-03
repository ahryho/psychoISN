library(SNPRelate)
library(data.table)

setwd("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/input/snps/")

# 1. Transform and save BED files using SNPRelate package

snp_bed_fn <- "bed/dex_geno_imputed.bed"
snp_bim_fn <- "bed/dex_geno_imputed.bim"
snp_fam_fn <- "bed/dex_geno_imputed.fam"

out_gds_fn <- "gds/dex_geno_imputed.gds"

snpgdsBED2GDS(bed.fn = snp_bed_fn, fam.fn = snp_fam_fn, bim.fn = snp_bim_fn, out_gds_fn)

# 2. Transform and save VCF file using SNPRelate package

snp_vcg_fn <- "dex_geno_imputed.vcf"

snpgdsVCF2GDS(vcf.fn = snp_vcg_fn, "gds/dex_geno_imputed_original.gds")