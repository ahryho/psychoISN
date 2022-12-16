#!/bin/bash

#SBATCH --job-name=mnda
#SBATCH --output=logs/mnda.out
#SBATCH --error=logs/mnda.err
#SBATCH --mem=500Gb
#SBATCH --part=pe
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=anastasiia_hry@psych.mpg.de

module load anaconda/anaconda3
module load R

eval "$(conda shell.bash hook)"

conda activate smccnet_psycho

rslt_dir=/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/results/5_fold_cv/all/mnda

input_graph_data=${rslt_dir}/global_network_mnda_input.csv
out_fn=${rslt_dir}/global_netqork_mnda_distances.rds

# Run MNDA

Rscript --vanilla "/home/ahryhorzhevska/kul/dex-stim-human-array-isns/code/05_mnda/02_mnda.R \
$input_graph_data $out_fn