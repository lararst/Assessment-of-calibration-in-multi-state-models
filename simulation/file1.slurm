#!/bin/bash
#SBATCH --job-name=file1
#SBATCH -o /output/dataset0
#SBATCH --export=ALL
#SBATCH --partition=cpu-zen4
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --mail-user=
#SBATCH --mail-type=END,FAIL

export R_LIBS="/R_pack"
module load gcc
module load  R/4.4.2-gfbf-2024a

Rscript /file1.R