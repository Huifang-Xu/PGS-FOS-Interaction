#!/bin/bash
#SBATCH --job-name=genoQC         # Job name
#SBATCH --partition=highmem_p               # Partition name (batch, highmem_p, or gpu_p)
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=20             # Number of CPU cores per task
#SBATCH --mem=900G                       # Memory per node (4GB); by default using M as unit
#SBATCH --time=23:00:00               # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=%x_%j.out              # Standard output log, e.g., testBowtie2_12345.out
#SBATCH --error=%x_%j.err               # Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=hx37930@uga.edu    # Where to send mail
#SBATCH --mail-type=END,FAIL # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --array=0-3

ml PLINK/1.9b_6.21-x86_64

pop=("EUR" "EAS" "CSA" "AFR")


p=${pop[${SLURM_ARRAY_TASK_ID}]}

imp_pass_File="ukb_mfi_chr${chr}_v3.INFO3.txt"
inDir="plink1.9_format"

#---------
#Set which
#steps run
#---------
step1=true
step2=true
step3=true
#---------

#######################################################
############# STEP 1. Genotype QC #####################
#######################################################
# Details:
# --geno remove variants with a missing genotype call rate > 0.05
# --mind 0.05 remove samples with a missing sample call rate > 0.05
# --mac 5 remove variants with a minor allele count < 5
# --hwe 1e-6 remove variants with a Hardy-Weinberg equilibrium P value < 1 × 10−6
#--keep Keep individuals that passed phenotype QC
# --exclude remove variants with an imputation quality < 0.3
if [ $step1 = true ]; then

echo "-=-=-=-=-=-=-=-STEP 1-=-=-=-=-=-=-=-\n\n"
# extract SNPs with an imputation quality INFO < 0.3

for p in ${pop[@]}
do
	outDir="bfile/${p}"
	keep_fam="phenotype/pgs/sampleID_${p}.txt"
	outputFile="${outDir}/ukb_chr${chr}.QCed"
	plink --bfile ${inDir}/ukb_chr${chr} \
	--geno 0.05 \
	--maf 0.001 \
	--hwe 1e-8 \
	--allow-no-sex \
	--exclude $imp_pass_File \
	--make-bed \
	--keep $keep_fam \
	--out $outputFile

done
fi

#######################################################
######## STEP 2. remove multiallelic variants #########
#######################################################
if [ $step2 = true ]; then

echo "-=-=-=-=-=-=-=-STEP 2-=-=-=-=-=-=-=-\n\n"
for p in ${pop[@]}
do
        outDir="bfile/${p}"
	cd ${outDir}
	# merge variants
	awk 'BEGIN{FS=OFS="\t"}{print $0,$1":"$4}' ukb_chr1.QCed.bim > ukb_allCHR.QCed.bim
	for c in {2..22};do awk 'BEGIN{FS=OFS="\t"}{print $0,$1":"$4}' ukb_chr${c}.QCed.bim >> ukb_allCHR.QCed.bim;done
	# check multiallelic variants
	ml R/4.3.1-foss-2022a
	Rscript scripts/sBayesRC/subset_biallelic.R ${outDir} 
	for c in {1..22}
	do
        	plink --bfile ${outDir}/ukb_chr${c}.QCed \
                	--extract ${outDir}/extract_biallelic.txt \
	                --make-bed \
	                --out ${outDir}/ukb_chr${c}.QCed.biallelic
	done
done
fi

#######################################################
########### STEP 3. merge all CHR files       #########
#######################################################
if [ $step3 = true ]; then

echo "-=-=-=-=-=-=-=-STEP 3-=-=-=-=-=-=-=-\n\n"
for p in ${pop[@]}
do
	outDir="bfile/${p}"
	plink --bfile ${outDir}/ukb_chr1.QCed.biallelic \
		--merge-list ${outDir}/bed_list.txt \
		--make-bed \
		--mind 0.05 \
		--out ${outDir}/ukb_allCHR_auto_QCed
done
fi
