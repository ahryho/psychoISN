#!/bin/bash

#SBATCH --job-name=isns_dist
#SBATCH --output=logs/veh_isns_dist_%A.out
#SBATCH --error=logs/veh_isns_dist_%A.err
#SBATCH --mem=400Gb
#SBATCH --part=pe
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=anastasiia_hry@psych.mpg.de

# if [ $# -ne 2 ]; then
#     echo "$0: usage: Not enough arguments
#           First argument: treatment (veh, dex, delta)"
#     exit 1
# fi

ulimit -s unlimited

module load R

dir_prefix=/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns

treatment=veh
chrom=all
cv_k=5
rslt_dir=${dir_prefix}/results/${cv_k}_fold_cv/$chrom/

Rscript --vanilla "/home/ahryhorzhevska/kul/dex-stim-human-array-isns/code/04_isns_analysis/02_distances_between_isns.R" \
$treatment $rslt_dir