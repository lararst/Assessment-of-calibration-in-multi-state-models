#!/bin/bash
#SBATCH --job-name=400mod_C1_s1  #C2/C3, 500/700
#SBATCH -o /output/400modC1_s1  #C2/C3, 500/700
#SBATCH --export=ALL
#SBATCH --partition=cpu-zen4
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --mail-user=
#SBATCH --mail-type=END,FAIL

export R_LIBS="/home/s4201302/R_pack"
module load gcc
module load  R/4.4.2-gfbf-2024a

Rscript /mod_C1_400_set1.R    #C2/C3, 500/700, set=1-50
