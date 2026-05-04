#!/bin/bash
#SBATCH --job-name=LDref_O3_Step1
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=80GB
#SBATCH --time=168:00:00

ml R/4.3.2-foss-2022b
export PATH=/home/gy71651/miniconda3/bin:$PATH

##############################################
# Variables: need to be fixed
ma_file="Omega_3_INT_COJO"                # GWAS summary data in COJO format
genotype="EUR/ukb_allCHR_auto_QCed"  # genotype prefix as LD refer EUR population
outDir="LD"             # Output folder that would be created aut$
threads=4                        # Number of CPU cores for eigen decomposition
#---usually don't need change bellow
genoCHR=""                       # If more than 1 genotype file, input range (e.g. "1-22") here.
refblock=""                      # Text file to define LD blocks, by default to use our GRCH37 coordination
tool="/home/gy71651/miniconda3/bin/gctb"                      # Command line to run gctb for generating the f$
start_idx=$1                     # Start index for this job
end_idx=$2                       # End index for this job
threads=4                        # Number of CPU cores for eigen decompos$

##############################################
# Code
# Step1: generate the LD block information and script
# Output $outDir/ldm.info, $outDir/ld.sh, $outDir/snplist/*.snplist
Rscript -e "SBayesRC::LDstep1(mafile='$ma_file', genoPrefix='$genotype', \
            outDir='$outDir', genoCHR='$genoCHR', blockRef='$refblock', log2file=TRUE)"

# Step 2
for idx in $(seq $start_idx $end_idx); do
    expected_file="$outDir/b${idx}.ldm.full.bin"
    if [[ ! -f $expected_file ]]; then
        echo "File $expected_file not found. Running LDstep2 for blockIndex $idx..."
        Rscript -e "SBayesRC::LDstep2(outDir='$outDir', blockIndex=$idx, log2file=TRUE)"
    else
        echo "File $expected_file already exists. Skipping blockIndex $idx."
    fi
done

# Step 3
for idx in $(seq $start_idx $end_idx); do
        Rscript -e "SBayesRC::LDstep3(outDir='$outDir', blockIndex=$idx, log2file=TRUE)"
done

# Step 4
Rscript -e "SBayesRC::LDstep4(outDir='$outDir', log2file=TRUE)"
