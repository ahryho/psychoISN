library(SNPRelate)
library(SmCCNet)
library(data.table)

snp_mpip <- fread("../../../dex-stim-human-array/data/integrative/matrixEQTL/snp_mtrx.csv")

setwd("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/input/snps")

# 1. Define the names 

snp_bed_fn <- "dex_geno_imputed.bed"
snp_bim_fn <- "dex_geno_imputed.bim"
snp_fam_fn <- "dex_geno_imputed.fam"

snp_vcg_fn <- "dex_geno_imputed.vcf"

# 2. Transform plink files using SNPRelate package

snpgdsBED2GDS(bed.fn = snp_bed_fn, fam.fn = snp_fam_fn, bim.fn = snp_bim_fn, "gds/dex_geno_imputed.gds")

# snpgdsVCF2GDS(vcf.fn = snp_vcg_fn, "gds/dex_geno_imputed_original.gds")

# 3. Explore transformed object

snps_gds <- snpgdsOpen("gds/dex_geno_imputed.gds")

# snpset <- snpgdsLDpruning(snps_gds, ld.threshold=0.2, maf = 0.05, slide.max.bp = 1000000L, slide.max.n = 100, method = "corr")

geno_obj   <- snpgdsGetGeno(snps_gds, with.id = T)

snpgdsClose(snps_gds)

snp_mtrx   <- geno_obj$genotype
colnames(snp_mtrx)   <- geno_obj$snp.id
rownames(snp_mtrx)   <- geno_obj$sample.id

head(snp_mtrx[, 1:10])
summary(snp_mtrx[, 1:10])

snp_mtrx <- snp_mtrx[rownames(snp_mtrx) %in% colnames(snp_mpip)[-1],]
dim(snp_mtrx)

snp_mtrx[is.na(snp_mtrx[, 2]), 2]
rs <- "rs3131971"
table(t(snp_mpip[SNP == rs, -1]))
table(snp_mtrx[, rs])

rs <- colnames(snp_mtrx)[7894]
table(t(snp_mpip[SNP == rs, -1]))
table(snp_mtrx[, rs])
