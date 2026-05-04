#! /usr/bin/Rscript
# R/4.3.1-foss-2022a
#This script is to munge GWAS summary statistics using MungeSumstat

Args <- commandArgs()

# Load required packages
#library(optparse)
library(MungeSumstats)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(data.table)
library(dplyr)

workDir <- Args[6]
inFile <- paste(workDir,Args[7],sep="/")
outFile <- paste(workDir,Args[8],sep="/")

df <- MungeSumstats::format_sumstats(path=inFile,
		      ref_genome="GRCh37",
#		      convert_small_p=TRUE, 
#		      allele_flip_check=TRUE,
#		      snp_ids_are_rs_ids=FALSE,
#		      INFO_filter=0.3,
		      return_data=TRUE,
		      nThread=2,
		      log_mungesumstats_msgs=TRUE, 
		      log_folder=workDir
		      )

#write result
write.table(df$sumstats,file=outFile,quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)

