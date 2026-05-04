#This folder contains scripts to calculate polygenic scores (PGS) of circulating omega-3 fatty acids


step0.extract_geno.sh: extraction and QC of UKB genotyping data 

MungSumstats.r: munge GWAS summary statistics using MungeSumstat 

step1.SBayesRC_munge_sumstats.nonUKB.sh: munge GWAS summary statistics from Karjalainen, M.K., et al., Nature, 2024 

step2.SBayesRC_PGS.nonEUR.sh: PGS calculation of UKB participants 

p_plus_t.clump_pgs.sh: script to perform LD clumping and PGS calculation
