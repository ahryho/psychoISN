library(SmCCNet)
library(gdsfmt)
library(SNPRelate)

# 1. Load DNAm beta mtrx

treatment <- "veh"
gds_fn   <- paste0("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/input/dnam/gds/methyl_beta_mtrx_corrected_for_cov",  
                   "_", treatment, ".gds") 

dnam_gds_file <- openfn.gds(gds_fn)

dnam_mtrx <- read.gdsn(index.gdsn(dnam_gds_file, "beta_mtrx"))
colnames(dnam_mtrx)   <- read.gdsn(index.gdsn(dnam_gds_file, "cpg_id"))
rownames(dnam_mtrx)   <- read.gdsn(index.gdsn(dnam_gds_file, "sample_id"))

closefn.gds(dnam_gds_file)

# 2. Load SNP data

snps_gds <- snpgdsOpen("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/input/snps/gds/dex_geno_imputed.gds")

# snpset <- snpgdsLDpruning(snps_gds, ld.threshold=0.2, maf = 0.05, slide.max.bp = 1000000L, slide.max.n = 100, method = "corr")

geno_obj  <- snpgdsGetGeno(snps_gds, with.id = T)
snp_chr   <- read.gdsn(index.gdsn(snps_gds, "snp.chromosome"))

snpgdsClose(snps_gds)

snp_mtrx   <- geno_obj$genotype
colnames(snp_mtrx)   <- geno_obj$snp.id
rownames(snp_mtrx)   <- geno_obj$sample.id

# 3. Load pheno trait