#!/bin/bash
#SBATCH --mem=100G
#SBATCH -J "HMM_prd"
#SBATCH --time=48:30:00
#SBATCH -o ../log/imp_%a_fs.out
#SBATCH -e ../log/imp_%a_fs.err
#SBATCH --partition=long
#SBATCH --array=1-3,5-6

ml R/4.1.0-foss-2021a
export R_MAX_VSIZE=128000000000

Rscript 04_PrdMisTP.R
