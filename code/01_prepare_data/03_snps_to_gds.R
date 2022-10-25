library(SNPRelate)
library(data.table)

setwd("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/input/snps/")

# 1. Transform and save BED files using SNPRelate package

snp_bed_fn <- "dex_geno_imputed.bed"
snp_bim_fn <- "dex_geno_imputed.bim"
snp_fam_fn <- "dex_geno_imputed.fam"

out_gds_fn <- "gds/dex_geno_imputed.gds"

snpgdsBED2GDS(bed.fn = snp_bed_fn, fam.fn = snp_fam_fn, bim.fn = snp_bim_fn, out_gds_fn)

# 2. Transform and save VCF file using SNPRelate package

snp_vcg_fn <- "dex_geno_imputed.vcf"

snpgdsVCF2GDS(vcf.fn = snp_vcg_fn, "gds/dex_geno_imputed_original.gds")

# 3. Transform and save Oxford .gen file using SNPRelate package

snp_gen_fn <- "plink2_gen/dex_geno_imputed_maf.gen"
sample_fn  <- "plink2_gen/dex_geno_imputed_maf.sample" 
out_gds_fn <- "gds/dex_geno_imputed_maf_from_gen.gds"

snpgdsGEN2GDS(gen.fn = snp_gen_fn, sample.fn = sample_fn, 
              out.fn = out_gds_fn, chr.code = 1, call.threshold = 0)

## The issue is that that SNPRelate package stores in "snp.id" the first column of GEN file which as a chromosome number
## and snp ids in a variabled called "snp.rs.id"
## In the next part, we'll re-write data so the var "snp.id" will contains SNP rs id.

snps_gds_file         <- openfn.gds(out_gds_fn, readonly = F)

delete.gdsn(index.gdsn(snps_gds_file, "snp.id"))
rename.gdsn(index.gdsn(snps_gds_file, "snp.rs.id"), newname = "snp.id")

closefn.gds(snps_gds_file)

snpgdsVCF2GDS(vcf.fn = snp_vcg_fn, "gds/dex_geno_imputed_original.gds")