#!/bin/bash
#
#SBATCH --job-name=smccnet_cv
#SBATCH --output=logs/smccnet_get_omic_modules_%A.out
#SBATCH --error=logs/smccnet_get_omic_modules_%A.err
#SBATCH --mem=900Gb    
#SBATCH --part=pe
#SBATCH --nodelist=pe7
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=anastasiia_hry@psych.mpg.de

# if [ $# -ne 2 ]; then
#     echo "$0: usage: Not enough arguments
#           First argument: treatment (veh, dex, delta)"
#     exit 1
# fi

module load R

treatment="veh"

Rscript --vanilla "/home/ahryhorzhevska/kul/dex-stim-human-array-isns/code/02_smccnet/02_smccnet_get_omics_modules.R" ${treatment}
