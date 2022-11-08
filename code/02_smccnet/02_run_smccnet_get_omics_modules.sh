#!/bin/bash

#SBATCH --job-name=smccnet
#SBATCH --output=logs/chromosomes/smccnet_%A_%a.out
#SBATCH --error=logs/chromosomes/smccnet_%A_%a.err
#SBATCH --mem-per-cpu=100Gb     # Each task uses max 100 Gb of memory
#SBATCH --array=1-21         # Submit 22 tasks. Run max 22 concurrently
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

treatment=veh
chrom=$(($SLURM_ARRAY_TASK_ID+1))
cv_k=5
cv_dir=${dir_prefix}/results/${cv_k}_fold_cv/chromosomes/${chrom}

dnam_gds_fn=${dir_prefix}/input/dnam/mad_filtered/gds/chromosomes/${treatment}/methyl_beta_mtrx_corrected_for_cov_${treatment}_chr${chrom}.gds
snps_gds_fn=${dir_prefix}/input/snps/ld_pruned/gds/chromosomes/dex_geno_chr${chrom}.gds

# Get omic modules

Rscript --vanilla "/home/ahryhorzhevska/kul/dex-stim-human-array-isns/code/02_smccnet/02_smccnet_get_omics_modules.R" \
$treatment $chrom $cv_k $cv_dir $dnam_gds_fn $snps_gds_fn