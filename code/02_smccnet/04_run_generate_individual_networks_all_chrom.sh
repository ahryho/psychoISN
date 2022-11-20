#!/bin/bash

#SBATCH --job-name=get_dex_individual_net_all_chr
#SBATCH --output=logs/get_dex_individual_net_all_chr_%A.out
#SBATCH --error=logs/get_dex_individual_net_all_chr_%A.err
#SBATCH --mem=800Gb
#SBATCH --part=pe
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=anastasiia_hry@psych.mpg.de

# if [ $# -ne 2 ]; then
#     echo "$0: usage: Not enough arguments
#           First argument: treatment (veh, dex, delta)"
#     exit 1
# fi

module load R

dir_prefix=/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns

treatment=dex
chrom=all
cv_k=5
cv_dir=${dir_prefix}/results/${cv_k}_fold_cv/$chrom/$treatment

dnam_gds_fn=${dir_prefix}/input/dnam/mad_filtered/gds/methyl_beta_mtrx_corrected_for_cov_mad80_filtered_${treatment}.gds
snps_gds_fn=${dir_prefix}/input/snps/ld_pruned/gds/dex_geno_imputed_maf_ld_pruned_from_gen.gds

echo Processing all chromosomes at once, i.e. no splitting

# Run SmCCNet CV

Rscript --vanilla "/home/ahryhorzhevska/kul/dex-stim-human-array-isns/code/02_smccnet/04_generate_individual_networks.R" \
$treatment $chrom $cv_dir $dnam_gds_fn $snps_gds_fn