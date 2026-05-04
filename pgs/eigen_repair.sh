#!/bin/bash
#SBATCH --job-name=LDref_O3_Step3
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=120GB
#SBATCH --time=168:00:00
#SBATCH --array=182,216,221,229,255,284,371,428,432,546

ml R/4.3.2-foss-2022b
export PATH=/home/gy71651/miniconda3/bin:$PATH

##############################################
# Variables
ma_file="Omega_3_INT_COJO"                # GWAS summary data
genotype="EUR/ukb_allCHR_auto_QCed"  # genotype file
outDir="LD"  # Output folder that stores the LD $
threads=4                        # Number of CPU cores for eigen decomposition
tool="/home/gy71651/miniconda3/bin/gctb"                      # Command line to run gctb

##############################################
Rscript eigen.R $SLURM_ARRAY_TASK_ID

