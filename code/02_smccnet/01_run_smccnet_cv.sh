#!/bin/bash
#
#SBATCH --job-name=smccnet_cv
#SBATCH --output=smccnet_cv_%A_%a.out
#SBATCH --error=smccnet_cv_%A_%a.err
#SBATCH --mem-per-cpu=2Gb     # Each task uses max 2Gb of memory
#SBATCH --part=pe
#SBATCH --exclude=hp01,hp02,pe6
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=anastasiia_hry@psych.mpg.de

# if [ $# -ne 2 ]; then
#     echo "$0: usage: Not enough arguments
#           First argument: treatment (veh, dex, delta)"
#     exit 1
# fi

module load R

treatment="veh"

Rscript --vanilla "/home/ahryhorzhevska/kul/dex-stim-human-array-isns/code/02_smccnet/01_smccnet_cv.R" ${treatment}
