#!bin/Rscript
library(dplyr)
library(data.table)

# accept parameters
Args <- commandArgs()

workDir = Args[6]
bim_all <- fread(paste0(workDir,"/ukb_allCHR.QCed.bim"), header=F)

dup_chrPos <- bim_all$V7[duplicated(bim_all$V7)]
dup_rsID <- bim_all$V2[duplicated(bim_all$V2)]

# remove duplicates
bim_bi <- subset(bim_all, ! (V2 %in% dup_rsID | V7 %in% dup_chrPos))

write.table(bim_bi[,1:6], paste0(workDir,"/ukb_allCHR.QCed.biallelic.bim"),row.names=F, col.names=F, quote=F,sep="\t")
write.table(bim_bi[,2], paste0(workDir,"/extract_biallelic.txt"),row.names=F, col.names=F, quote=F)

