library(SNPRelate)
library(data.table)
library(arrow)

setwd("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/input/snps/")

# 1. Read arguments

args       <- commandArgs(trailingOnly=TRUE)
snp_gen_fn <- args[1]
sample_fn  <- args[2] 
out_gds_fn <- args[3] # "gds/dex_geno_imputed_maf_from_gen.gds"

# 2. Transform and save Oxford .gen file using SNPRelate package

snpgdsGEN2GDS(gen.fn = snp_gen_fn, sample.fn = sample_fn, 
              out.fn = out_gds_fn, chr.code = 1, call.threshold = 0)

## 2.1. Change SNP id names

## The issue is that that SNPRelate package stores in "snp.id" the first column of GEN file which as a chromosome number
## and snp ids in a variabled called "snp.rs.id"
## In the next part, we'll re-write data so the var "snp.id" will contains SNP rs id.

snps_gds_file         <- openfn.gds(out_gds_fn, readonly = F)

delete.gdsn(index.gdsn(snps_gds_file, "snp.id"))
rename.gdsn(index.gdsn(snps_gds_file, "snp.rs.id"), newname = "snp.id")

## 2.2. Add chrom and pos into gds file

snp_loc_fn <- "snp_locations.csv"
snp_loc    <- read_delim_arrow(snp_loc_fn, delim = ";")

snp_ids    <- read.gdsn(index.gdsn(snps_gds_file, "snp.id"))

order_idx  <- match(snp_ids, snp_loc$SNP)
snp_loc    <- snp_loc[order_idx, ]
all(snp_loc$SNP == snp_ids)

delete.gdsn(index.gdsn(snps_gds_file, "snp.chromosome"))
add.gdsn(snps_gds_file, name = "snp.chromosome", val = snp_loc$chr)

print("Check GDS file: ", quote = F)
snps_gds_file

closefn.gds(snps_gds_file)
