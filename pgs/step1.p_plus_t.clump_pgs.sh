#!/bin/bash
#SBATCH --job-name=pgs_nonUKB         # Job name
#SBATCH --partition=batch               # Partition name (batch, highmem_p, or gpu_p)
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=150G                       # Memory per node (4GB); by default using M as unit
#SBATCH --time=1-23:00:00               # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=%x_%j.out              # Standard output log, e.g., testBowtie2_12345.out
#SBATCH --error=%x_%j.err               # Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=hx37930@uga.edu    # Where to send mail
#SBATCH --mail-type=END,FAIL # Mail events (BEGIN, END, FAIL, ALL)

ml PLINK/1.9b_6.21-x86_64

trait_ids=("DHA" "DHA_pct" "Omega_3" "Omega_3_pct")

inDir="gwasSumstats_nonUKB_munged"
outDir="10.P_plus_T"
bfile="/bfile/${p}/ukb_allCHR_auto_QCed"
glist_hg19="glist-hg19"

# step 1: Flip A2 as effect allele, flip beta and freq to A2 as effect allele
for t in {0..3}
do
        trait=${trait_ids[$t]}
	inFile=${inDir}/${trait}.A1effect.munged.txt
	clumpedFile=${outDir}/clumped/${trait}
	clumped_sumstat=${outDir}/clumped/${trait}.clumped.sumstat.txt
	pgsFile=${outDir}/pgs/${trait}
	# flip A2 as A1 effect allele
	awk 'BEGIN{FS=OFS="\t"}NR==1{print "SNP","CHR","BP","A1","A2","BETA","SE","P","FRQ","DIRECTION","HETISQT","HETCHISQ","HETDF","HETPVAL"}NR>1{print $1,$2,$3,$5,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14}' ${inDir}/${trait}.A2effect.munged.txt > ${inDir}/${trait}.A1effect.munged.txt
	# clumping
	plink --bfile ${bfile} --clump ${inFile} --clump-range ${glist_hg19} --out ${clumpedFile} --allow-extra-chr --clump-p1 5e-8 --clump-p2 5e-8 --clump-r2 0.1 --clump-kb 250 --clump-field P --clump-snp-field SNP --threads 20
	# Extract clumped SNPs
	awk 'NR>1{print $3}' ${clumpedFile}.clumped |grep "rs"  > ${clumpedFile}.clumped.snp
	head -n 1 ${inFile} > ${clumped_sumstat}
	for i in `cat ${clumpedFile}.clumped.snp`; do echo "awk 'BEGIN{FS=OFS=\"\\t\"}NR>1{if(\$1==\"${i}\")print}' ${inFile} >> ${clumped_sumstat}" >> tmp.sh;done
	sh tmp.sh; rm tmp.sh
	# calculate PGS
	plink --bfile ${bfile} --score ${clumped_sumstat} 1 4 6 header --out ${pgsFile}
done
