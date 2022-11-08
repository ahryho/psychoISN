# LD pruning on MAFed genotypes

cd /binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/input/snps/ld_pruned

# Run LD pruning

plink --bfile ../unpruned/bed/dex_geno_imputed_maf --exclude range mhc_regions.txt \
      --indep-pairwise 100 50  0.2 --out dex_geno_imputed_maf_ld
      
mv dex_geno_imputed_maf_ld.log logs

# Extract pruned SNPs

plink --bfile ../unpruned/bed/dex_geno_imputed_maf --extract dex_geno_imputed_maf_ld.prune.in \
      --make-bed --out bed/dex_geno_imputed_maf_ld_pruned

mv bed/dex_geno_imputed_maf_ld_pruned.log logs

# Generate oxford genome file

mkdir -p gen
plink --bfile bed/dex_geno_imputed_maf_ld_pruned --recode oxford --out gen/dex_geno_imputed_maf_ld_pruned

mv gen/dex_geno_imputed_maf_ld_pruned.log logs/dex_geno_imputed_maf_ld_pruned_gen.log

# Generate GDS object from Oxford .gen

snp_gen_fn="ld_pruned/gen/dex_geno_imputed_maf_ld_pruned.gen"
sample_gen_fn="ld_pruned/gen/dex_geno_imputed_maf_ld_pruned.sample"
out_gds_fn="ld_pruned/gds/dex_geno_imputed_maf_ld_pruned_from_gen.gds"

Rscript --vanilla ~/kul/dex-stim-human-array-isns/code/00_functions/gen2gds.R $snp_gen_fn $sample_gen_fn $out_gds_fn
