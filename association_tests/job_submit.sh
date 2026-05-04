#!/bin/bash
#SBATCH --job-name=asso         # Job name
#SBATCH --partition=batch               # Partition name (batch, highmem_p, or gpu_p)
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=10             # Number of CPU cores per task
#SBATCH --mem=120G                       # Memory per node (4GB); by default using M as unit
#SBATCH --time=3-23:00:00               # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=%x_%j.out              # Standard output log, e.g., testBowtie2_12345.out
#SBATCH --error=%x_%j.err               # Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=hx37930@uga.edu    # Where to send mail
#SBATCH --mail-type=END,FAIL # Mail events (BEGIN, END, FAIL, ALL)

ml R/4.3.1-foss-2022a

p="EUR"
pgs_suffix=".profile"
phenoFile="pheno_Sep2025.clean.fullPop.QCed.FOS.Oily.NMR.imputed.tsv"
inDir="10.P_plus_T/pgs/${p}"
outDir="10.P_plus_T/asso_result/${p}"
echo "
===================================================
=================== ${p}: primary EUR =============
==================================================="
Rscript asso_pgs_pheno.primary.r ${p} ${pgs_suffix} ${phenoFile} ${inDir} ${outDir}

pops=("EUR" "AFR" "EAS" "CSA")
phenoFile="pheno_Nov2025.clean.fullPop.QCed.exEURr1r2.imputed.FOS.Oily.tsv"
outDir="10.P_plus_T/asso_result_rep/${p}"
for p in ${pops[@]}
do
echo "
===================================================
== ${p}: replication (EUR) and primary (non-EUR) ==
==================================================="	
Rscript asso_pgs_pheno.primary.r ${p} ${pgs_suffix} ${phenoFile} ${inDir} ${outDir}
done

p="EUR"
phenoFile="pheno_Sep2025.clean.fullPop.QCed.imputed.24hDiet.all.tsv"
outDir="10.P_plus_T/asso_result_24hDiet/${p}"
echo "
===================================================
=================== ${p}: 24h diet ================
==================================================="
        Rscript asso_pgs_pheno.diet24h.r ${p} ${pgs_suffix} ${phenoFile} ${inDir} ${outDir}

