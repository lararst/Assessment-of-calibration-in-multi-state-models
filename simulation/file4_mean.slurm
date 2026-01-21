#!/bin/bash
#SBATCH --job-name=mean_C1_400        #C2/C3, 500/700
#SBATCH -o /output/C1_400    #C2/C3, 500/700
#SBATCH --export=ALL
#SBATCH --partition=cpu-zen4
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-user=
#SBATCH --mail-type=END,FAIL

export R_LIBS="/R_pack"
module load gcc
module load  R/4.4.2-gfbf-2024a

Rscript /mean_C1_400.R        #C2/C3, 500/700