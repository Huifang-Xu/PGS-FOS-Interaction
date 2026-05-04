'%ni%' <- Negate('%in%')
library(plyr)
library(dplyr)
library(tidyverse)
library(readr)
library(data.table)


#################################################################################
###################### Read data: full population ###############################
#################################################################################
# read in phenotype file, including 251 metabolites from release 3
pheno <- fread("../NMR_release3/NMR_metabolomics_rawdata_50w.csv",header=T, sep=",")
pheno <- as_tibble(pheno)
pheno$eid <- as.character(pheno$eid)
old_colnames <- names(pheno)
new_colnames <- paste0("f.",old_colnames)
new_colnames <- gsub(pattern = "-", replacement = ".",new_colnames)
names(pheno) <- new_colnames

# read in phenotype file that is analyzed in the discovery analysis.
pheno_pri <- fread("pheno_Sep2025.clean.fullPop.QCed.FOS.Oily.NMR.imputed.tsv")
pheno_pri <- as_tibble(pheno_pri)
pheno_pri$f.eid <- as.character(pheno_pri$f.eid)
# exclude EUR participants who are already analyzed
exclude_ids <- pheno_pri$f.eid[pheno_pri$pop == "EUR"]

# read in QCed covariate file
covar <- fread("pheno_Sep2025.clean.fullPop.QCed.tsv")
covar <- as_tibble(covar)
covar$f.eid <- as.character(covar$f.eid)
covar <- covar[,c(1:313,347:389)]

# read in imputed covariate file
imputed_data <- readRDS("pheno_Sep2025.clean.fullPop.QCed.IMPUTE_m5_maxit20.rds")
complete_data <- complete(imputed_data,1)
names(complete_data) <- c('f.eid', 'sex', 'age', 'BMI_imp', 'TSI_imp', 'smoking_imp', 'alcohol_imp', 'physical_imp','Statin','pop')
complete_data$f.eid <- as.character(complete_data$f.eid)
# merge with covariate data
covar_imputed <- full_join(covar,complete_data[,c('f.eid','BMI_imp', 'TSI_imp', 'smoking_imp', 'alcohol_imp', 'physical_imp')],by='f.eid')

# merge data
pheno_merge <- right_join(pheno,covar_imputed,by="f.eid")

# exclude EUR participants who are already analyzed (NMR metabolomics release 1 & 2)
pheno_rep <- subset(pheno_merge, ! f.eid %in% exclude_ids)

write.table(pheno_rep, "replication_NMR_release3/clean.fullPop.QCed.exEURr1r2.imputed.tsv",sep="\t",col.names=T,row.names=F,quote=F)

# remove participants without Omega-3 NMR measurement
pheno_rep <- pheno_rep[!is.na(pheno_rep$f.23444.0.0),]

# remove individuals without fish oil intake data
table(pheno_rep$f.6179.0.0)
#    -7     -3      1      2      3      4      5      6
#113585   1171  62453  14321   4937   1636   2401    495
pheno_rmFOS <- pheno_rep[!is.na(pheno_rep$f.6179.0.0),]
pheno_rmFOS <- pheno_rep[pheno_rep$f.6179.0.0 != -3,]
pheno_rmFOS$FishOil <- 0
pheno_rmFOS$FishOil[pheno_rmFOS$f.6179.0.0 ==1] <- 1

#################################################################################
# preprocessing oily fish intake: stratified oily fish intake into two groups: low intake group: never, less than once a week. High intake group: once a week, 2–4 times a week, 5–6 times a week, and once or more daily
#################################################################################
table(pheno_rmFOS$f.1329.0.0)
#   -3    -1     0     1     2     3     4     5
#  101  1224 22405 66243 74399 33382  1500   574
# stratified oily fish intake into two groups: low vs high intake.
pheno_rmFOS$OilyFish <- NA
pheno_rmFOS$OilyFish[pheno_rmFOS$f.1329.0.0 >= 2] <- 1
pheno_rmFOS$OilyFish[pheno_rmFOS$f.1329.0.0 == 0 | pheno_rmFOS$f.1329.0.0 ==1] <-0

write.table(pheno_rmFOS, "replication_NMR_release3/pheno_Sep2025.clean.fullPop.QCed.exEURr1r2.imputed.FOS.Oily.tsv",sep="\t",col.names=T,row.names=F,quote=F)

#################################################################################
# preprocessing 24h dietary recall: consider fish oil user as long as they answer yes in one of five instances. Prerequisite: participants must answer data field f.104670.X.X at each instance
#################################################################################
pheno_rep <- pheno_merge[pheno_merge$pop == "EUR",]
result <- pheno_rep[,'f.eid']
instances <- seq(0,4,1)

# precess each instance
for (i in instances) {
# filter participants without info on f.104670
colName <- paste0('f.104670.',i, '.0')
eligible_ids <- pheno_rep$f.eid[!is.na(pheno_rep[[colName]])]

# extract FOS of each instance
datafield1 <- paste0("f.20084.",i,".")
datafield2 <- paste0("f.20084.",i,".0")
tmp <- pheno_rep %>% select(f.eid,starts_with(datafield1))
response_ids <- tmp$f.eid[!is.na(tmp[[datafield2]])]

# recode FOS as 1 (Fish oil user) or 0 (non fish oil user). FOS code: 472
colN <- ncol(tmp)
tmp[,2:colN][tmp[,2:colN] != 472] <- 0
tmp[,2:colN][tmp[,2:colN] == 472] <- 1
# get FOS stat
new_colName <- paste0("FOS_24h_",i)
tmp[[new_colName]] <- rowSums(tmp[,2:colN],na.rm = TRUE)
# if participant missing information on f.104670, then change the value to NA
tmp[[new_colName]][! tmp$f.eid %in% eligible_ids] <- NA
# if participant have information on f.104670 but answer no to fish oil (!= 472), they are treated as non-user
tmp[[new_colName]][tmp[[new_colName]] !=1 & (tmp$f.eid %in% eligible_ids)] <- 0

# merge five instances
result <- left_join(result,tmp[,c('f.eid',new_colName)],by='f.eid')
}

# sum five instances
result$FOS_24h <- rowSums(result[,2:6],na.rm = TRUE)

# obtain info on eligible participants who fill in 24h dietary recall questionnair at least one
eligible_ids_24h <- unique(c(result$f.eid[!is.na(result$FOS_24h_0)], result$f.eid[!is.na(result$FOS_24h_1)],result$f.eid[!is.na(result$FOS_24h_2)],result$f.eid[!is.na(result$FOS_24h_3)],result$f.eid[!is.na(result$FOS_24h_4)]))

# recode FOS as 1 (Fish oil user) or 0 (non fish oil user)
result$FOS_24h[result$FOS_24h > 1] <- 1
result$FOS_24h[result$FOS_24h == 0 & ! (result$f.eid %in% eligible_ids_24h)] <- NA

# merge with other phenotype
pheno_24hDiet <- left_join(pheno_rep,result,by='f.eid')

# remove participants without information on 24h dietary recall
pheno_24hDiet_clean <- pheno_24hDiet[!is.na(pheno_24hDiet$FOS_24h),]

write.table(pheno_24hDiet_clean, "replication_NMR_release3/pheno_Sep2025.clean.EUR.QCed.allPhases.imputed.24hDiet.tsv",sep="\t",col.names=T,row.names=F,quote=F)

