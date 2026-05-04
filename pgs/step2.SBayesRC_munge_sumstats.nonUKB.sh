#!/bin/bash
#SBATCH --job-name=munge_nonUKB         # Job name
#SBATCH --partition=batch               # Partition name (batch, highmem_p, or gpu_p)
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=150G                       # Memory per node (4GB); by default using M as unit
#SBATCH --time=1-23:00:00               # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=%x_%j.out              # Standard output log, e.g., testBowtie2_12345.out
#SBATCH --error=%x_%j.err               # Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=hx37930@uga.edu    # Where to send mail
#SBATCH --mail-type=END,FAIL # Mail events (BEGIN, END, FAIL, ALL)


ml SNPlocs.Hsapiens.dbSNP155.GRCh37/0.99.24-foss-2024a-R-4.4.2

GSCT_ids=("GCST90301955" "GCST90301956" "GCST90301959" "GCST90301960")
trait_ids=("DHA" "DHA_pct" "Omega_3" "Omega_3_pct")

inDir="gwasSumstats_nonUKB"
outDir="gwasSumstats_nonUKB_munged/noFilter"
cojoDir="gwasSumstats_nonUKB_COJO"

# step 1: Flip A2 as effect allele, flip beta and freq to A2 as effect allele
for t in {0..3}
do
	GSCT=${GSCT_ids[$t]}
	trait=${trait_ids[$t]}
	# step 1: Flip A2 as effect allele, flip beta and freq to A2 as effect allele
	awk 'BEGIN{FS=OFS="\t"}NR==1{print "CHR","BP","A1","A2","BETA","SE","P","FRQ","direction","hetisq","hetchisq","hetdf","hetpval"}NR>1{print $1,$2,$4,$3,$5,$6,$8,$7,$10,$11,$12,$13,$14}' ${inDir}/${GSCT}.tsv > ${outDir}/${trait}.A2effect.txt
	# step 2: run MungSumstats: A2 is effect allele, beta and freq are refer to A2
	Rscript MungSumstats.r ${outDir} ${trait}.A2effect.txt ${trait}.A2effect.munged.txt
	# step 3: flip effect allele of A2 back to A1 so that all GWAS sumstat are treating the same alleles
	awk 'BEGIN{FS=OFS="\t"}NR==1{print "SNP","CHR","BP","A1","A2","BETA","SE","P","FRQ","DIRECTION","HETISQT","HETCHISQ","HETDF","HETPVAL"}NR>1{print $1,$2,$3,$5,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14}' ${outDir}/${trait}.A2effect.munged.txt > ${outDir}/${trait}.A1effect.munged.txt
	# step4: convert to COJO format
	awk 'BEGIN{FS=OFS="\t"}{print $1,$4,$5,$9,$6,$7,$8,"136016"}' ${outDir}/${trait}.A1effect.munged.txt > ${cojoDir}/${trait}.A1effect.munged.COJO.txt
done
