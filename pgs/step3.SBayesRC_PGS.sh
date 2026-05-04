#!/bin/bash
#SBATCH --job-name=SBayesRC
#SBATCH --partition=iob_batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=500GB
#SBATCH --time=5-23:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=hx37930@uga.edu
#SBATCH --mail-type=END,FAIL

#ml R/4.3.2-foss-2022b
ml SBayesRC/0.2.6-foss-2024a-R-4.4.2
ml PLINK/2.0.0-a.6.20-gfbf-2024a

trait_ids=("DHA" "DHA_pct" "Omega_3" "Omega_3_pct")

pops=("EUR" "CSA" "AFR" "EAS")

ld_folder="09.SBayesRC/LD"
annot="scripts/annot_baseline2.2.txt"
outDir="09.SBayesRC/result_nonUKBgwas"
threads=20

export OMP_NUM_THREADS=$threads #
export PATH=/apps/eb/PLINK/2.0.0-a.6.20-gfbf-2024a/bin/plink2:$PATH

for t in {0..14}
do
trait=${trait_ids[$t]}
for p in ${pops[@]}
do
ma_file="09.SBayesRC/gwasSumstats_nonUKB_COJO/${trait}.A1effect.munged.COJO.txt"
out_prefix="09.SBayesRC/result_nonUKBgwas/${trait}/${trait}_QCed"

mkdir -p ${outDir}/${trait}

#---------
#Set which
#steps run
#---------
step1=true
step2=true
step3=true
step4=true
#---------

###############################################################################
########## STEP 1. QC  the  summary statistics with tidy function #############
###############################################################################
if [ $step1 = true ]; then
echo "-=-=-=-=-=-=-=-STEP 1-=-=-=-=-=-=-=-\n\n"
Rscript -e "SBayesRC::tidy(mafile='$ma_file', LDdir='$ld_folder', output='${out_prefix}_tidy.ma', freq_thresh=0.2,
	rate2pq=0.5, log2file=TRUE)"
fi

###############################################################################
########## STEP 2. Impute the ma_file to match the LD #########################
###############################################################################
if [ $step2 = true ]; then
echo "-=-=-=-=-=-=-=-STEP 2-=-=-=-=-=-=-=-\n\n"
Rscript -e "SBayesRC::impute(mafile='${out_prefix}_tidy.ma', LDdir='$ld_folder', 
	output='${out_prefix}_imp.ma', log2file=TRUE)"
fi

###############################################################################
########## STEP 3. Run main function sbayesrc #################################
###############################################################################
if [ $step3 = true ]; then
echo "-=-=-=-=-=-=-=-STEP 3-=-=-=-=-=-=-=-\n\n"
Rscript -e "SBayesRC::sbayesrc(mafile='${out_prefix}_imp.ma', LDdir='$ld_folder', 
outPrefix='${out_prefix}_sbrc', annot='$annot', log2file=TRUE)"
fi


###############################################################################
########## STEP 4. Generate PRS with sbayes output ############################
###############################################################################
if [ $step4 = true ]; then
echo "-=-=-=-=-=-=-=-STEP 4-=-=-=-=-=-=-=-\n\n"
genoPrefix="/bfile/${p}/ukb_allCHR_auto_QCed"
output="${out_prefix}_SBayesRC_PGS_${p}"
Rscript -e "SBayesRC::prs(weight='${out_prefix}_sbrc.txt', genoPrefix='$genoPrefix', out='$output')"
fi

done
done
