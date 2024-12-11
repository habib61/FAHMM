#!/bin/bash
#SBATCH --mem=3500
#SBATCH --time=3:30:00
#SBATCH --output="../log/fahmm_%a.out"
#SBATCH --error="../log/fahmm_%a.err"
#SBATCH --partition=long,main
#SBATCH --array=1-12

ml R/4.1.0-foss-2021a

Rscript 06_fahmm.R 
