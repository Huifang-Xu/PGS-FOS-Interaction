# Summary

This folder contains scripts for calculating polygenic scores (PGS) for circulating omega-3 fatty acids.

step0.extract_geno.sh: extraction and QC of UKB genotyping data \
step1.p_plus_t.clump_pgs.sh: script to perform LD clumping and PGS calculation using P+T approach

step1.SBayesRC_LDref.sh: Calculate in-sample LD reference from 239,268 UKB participants of EUR ancestry.

MungSumstats.r: munge GWAS summary statistics using MungeSumstat 

step2.SBayesRC_munge_sumstats.nonUKB.sh: munge GWAS summary statistics from Karjalainen, M.K., et al., Nature, 2024 

step3.SBayesRC_PGS.nonEUR.sh: PGS calculation of UKB participants 


