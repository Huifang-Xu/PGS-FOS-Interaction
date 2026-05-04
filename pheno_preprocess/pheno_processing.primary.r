'%ni%' <- Negate('%in%')
library(plyr)
library(dplyr)
library(tidyverse)
library(readr)
library(data.table)
library(rsq)

#################################################################################
###################### Read data: full population ###############################
#################################################################################
# read in relatedness file
rel <- fread("ukb_rel_a48818_s488117.dat",header=T,sep=" ")
rel <- as_tibble(rel)
rel$ID1 <- as.character(rel$ID1)
rel$ID2 <- as.character(rel$ID2)

# read in PCA covariates file
pc <- fread("all_pops_non_eur_pruned_within_pop_pc_covs.tsv",header=T,sep="\t")
pc <- as_tibble(pc)
pc$s <- as.character(pc$s)

# read in sample ID file
sampleID <- fread("ukb48818bridge31063.txt",header=F,sep=" ")
sampleID <- as_tibble(sampleID)
colnames(sampleID) <- c("f.eid","s")
sampleID$f.eid <- as.character(sampleID$f.eid)
sampleID$s <- as.character(sampleID$s)

# read in withdrwan samples
withdrawn <- fread("withdrawn_w48818_20250818.csv",header=F)
withdrawn <- as_tibble(withdrawn)
colnames(withdrawn) <- "f.eid"
withdrawn$f.eid <- as.character(withdrawn$f.eid)

# read in phenotype file, including 249 metabolites and covariates
pheno <- fread("pheno_Sept2025.txt",header=T, sep="\t")
pheno <- as_tibble(pheno)
pheno$f.eid <- as.character(pheno$f.eid)

# read in Townsend index file f.189.0.0
townsend <- fread("townsend.txt", header=TRUE, sep="\t")
townsend <- as_tibble(townsend)
townsend$f.eid <- as.character(townsend$f.eid)


#### function of combining multiple columns into one
manyColsToDummy<-function(search_terms, search_columns,output_table){
        #initialize output table
        temp_table<-data.frame(matrix(ncol=length(search_terms), nrow= nrow(search_columns)))
        colnames(temp_table)<-search_terms

        #make table
        for (i in 1:length(search_terms)){
                vec<-rowSums(sapply(search_columns,function(x) grepl(search_terms[i], x, ignore.case = TRUE)))>0
                temp_table[,i]<-vec
        }
        temp_table<-sapply(temp_table, as.integer, as.logical)
        temp_table<-as.data.frame(temp_table)
        assign(x = output_table, value = temp_table, envir = globalenv())
}


#################################################################################
#################### phenotype pre-processing: full population ##################
#################################################################################
# 1. merge datasets with Townsend index (TSI)
pheno <- left_join(pheno, townsend, by = 'f.eid')

# 2. genotyping array
pheno$f.22000.0.0=ifelse(pheno$f.22000.0.0>0,1,0)

# 3. Assessment_centres: code 10003 is missed in AFR participants,code 11023 is missed in CSA participants
pheno$f.54.0.0 <- mapvalues(as.character(pheno$f.54.0.0), c(11012, 11021,11011,11008,11003,11020,11005,11004,11018,11010,11016,11001,11017,11009,11013,11002,11007,11014,11006,11022,11023,10003), c("a11012","a11021","a11011","a11008","a11003","a11020",'a11005',"a11004","a11018","a11010","a11016","a11001","a11017","a11009","a11013","a11002","a11007","a11014","a11006","a11022","a11023","a10003"))
pheno$Assessment_centres_10003=ifelse(pheno$f.54.0.0=='a10003',1,0)
pheno$Assessment_centres_11001=ifelse(pheno$f.54.0.0=='a11001',1,0)
pheno$Assessment_centres_11002=ifelse(pheno$f.54.0.0=='a11002',1,0)
#non-England
pheno$Assessment_centres_11003=ifelse(pheno$f.54.0.0=='a11003',1,0)
#non-England
pheno$Assessment_centres_11004=ifelse(pheno$f.54.0.0=='a11004',1,0)
#non-England
pheno$Assessment_centres_11005=ifelse(pheno$f.54.0.0=='a11005',1,0)
pheno$Assessment_centres_11006=ifelse(pheno$f.54.0.0=='a11006',1,0)
pheno$Assessment_centres_11007=ifelse(pheno$f.54.0.0=='a11007',1,0)
pheno$Assessment_centres_11008=ifelse(pheno$f.54.0.0=='a11008',1,0)
pheno$Assessment_centres_11009=ifelse(pheno$f.54.0.0=='a11009',1,0)
pheno$Assessment_centres_11010=ifelse(pheno$f.54.0.0=='a11010',1,0)
pheno$Assessment_centres_11011=ifelse(pheno$f.54.0.0=='a11011',1,0)
pheno$Assessment_centres_11012=ifelse(pheno$f.54.0.0=='a11012',1,0)
pheno$Assessment_centres_11013=ifelse(pheno$f.54.0.0=='a11013',1,0)
pheno$Assessment_centres_11014=ifelse(pheno$f.54.0.0=='a11014',1,0)
pheno$Assessment_centres_11016=ifelse(pheno$f.54.0.0=='a11016',1,0)
pheno$Assessment_centres_11017=ifelse(pheno$f.54.0.0=='a11017',1,0)
pheno$Assessment_centres_11018=ifelse(pheno$f.54.0.0=='a11018',1,0)
pheno$Assessment_centres_11020=ifelse(pheno$f.54.0.0=='a11020',1,0)
pheno$Assessment_centres_11021=ifelse(pheno$f.54.0.0=='a11021',1,0)
#non-England
pheno$Assessment_centres_11022=ifelse(pheno$f.54.0.0=='a11022',1,0)
#non-England
pheno$Assessment_centres_11023=ifelse(pheno$f.54.0.0=='a11023',1,0)

# 4. merge sampleID; match s (PCA) with eid (phenotype)
sample_pca <- unique(subset(pc,select=c(s, pop, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)))
sample_pca <- sample_pca %>%  distinct(s, pop, .keep_all = TRUE)

sampleID_merge <- sample_pca %>% inner_join(sampleID,by="s")

pheno_mergeID <- pheno %>% inner_join(sampleID_merge,by="f.eid")

# 5. statin use
statin_cols<-c(sprintf("f.20003.0.%s", 0:47))
statin_codes<-c(1141146234,1141192414,1140910632,1140888594,1140864592,1141146138,1140861970,1140888648,1141192410,1141188146,1140861958,1140881748,1141200040)
# extract columns of statin medications
dt_statin <- pheno_mergeID %>% select(statin_cols)
dt_statin_ID <- pheno_mergeID %>% select(f.eid,statin_cols)
# data transformation of statin
manyColsToDummy(statin_codes, dt_statin, "bd_temp_statin")
bd_temp_statin$Statin <- rowSums(bd_temp_statin)
bd_temp_statin$Statin[bd_temp_statin$Statin>1] <- 1
bd_temp_statin$f.eid <- dt_statin_ID$f.eid
bd_temp_statin <- bd_temp_statin %>% select(f.eid, Statin)
# merge datasets
pheno_addStatin <- pheno_mergeID %>% left_join(bd_temp_statin, by="f.eid")

# add sex_by_age column
pheno_addStatin$sex_by_age <- pheno_addStatin$f.22001.0.0 * pheno_addStatin$f.21003.0.0

#write.table(pheno_mergeID, "pheno_249metabolites.clean.fullPop.tsv",sep="\t",col.names=T,row.names=F,quote=F)

##################################################################################
# phenotype QC (full population): 1) remove individuals with withdrawn info; 2) remove individuals with 3rd-degree relatives or closer; 3) mismatched information between phenotypic and genetix sex; 4) Sex chromosomes aneuploidy; 5) outlier for het and missing genotype rate
##################################################################################
# 1. remove samples with withdrawn info (full population)
pheno_rmWithdrwan <- subset(pheno_addStatin, ! f.eid %in% withdrawn$f.eid)
  
# 2. remove individuals with high degree of genetic kinship
pheno_rmKinship <- pheno_rmWithdrwan[pheno_rmWithdrwan$f.22021.0.0 != 10,]

# 3. mismatched information between phenotypic (f.31.0.0) and genetic sex (f.22001.0.0)
pheno_rmMismatchedSex <-  pheno_rmKinship[pheno_rmKinship$f.31.0.0 == pheno_rmKinship$f.22001.0.0,]

# 4. Sex chromosomes aneuploidy
pheno_rmSexAneuploidy <- pheno_rmMismatchedSex %>%filter(is.na(f.22019.0.0)) 

# 5. outlier for het and missing genotype rate
pheno_rmHetMiss <- pheno_rmSexAneuploidy %>%filter(is.na(f.22027.0.0))

write.table(pheno_rmHetMiss, "pheno_Sep2025.clean.fullPop.QCed.tsv",sep="\t",col.names=T,row.names=F,quote=F)

# read in imputed data
imputed_data <- readRDS("pheno_Sep2025.clean.fullPop.QCed.IMPUTE_m5_maxit20.rds")
complete_data <- complete(imputed_data,1)
names(complete_data) <- c('f.eid', 'sex', 'age', 'BMI_imp', 'TSI_imp', 'smoking_imp', 'alcohol_imp', 'physical_imp','Statin','pop')
# merge with raw data
pheno_rmHetMiss_imputed<- left_join(pheno_rmHetMiss, complete_data[,c('f.eid','BMI_imp', 'TSI_imp', 'smoking_imp', 'alcohol_imp', 'physical_imp')], by='f.eid')

#################################################################################
# remove individuals without fish oil intake data
#################################################################################
table(pheno_rmHetMiss$f.6179.0.0)
#    -7     -3      1      2      3      4      5      6
#252317   1203 140308  31726  10306   3597   4755   1161
pheno_rmFOS <- pheno_rmHetMiss[pheno_rmHetMiss$f.6179.0.0 != -3,]
# remove individuals with missing information on f.6179.0.0
#pheno_rmFOS <- pheno_rmHetMiss[!is.na(pheno_rmHetMiss$f.6179.0.0),]
pheno_rmFOS$FishOil <- 0
pheno_rmFOS$FishOil[pheno_rmFOS$f.6179.0.0 ==1] <- 1

write.table(pheno_rmFOS, "pheno_Sep2025.clean.fullPop.QCed.FOS.tsv",sep="\t",col.names=T,row.names=F,quote=F)

# keep participants with NMR measurements, use omega-3 as major outcome
pheno_rmFOS_NMR <- pheno_rmFOS[!is.na(pheno_rmFOS$f.23444.0.0),]

write.table(pheno_rmFOS_NMR, "pheno_Sep2025.clean.fullPop.QCed.FOS.NMR.tsv",sep="\t",col.names=T,row.names=F,quote=F)

# imputed data
pheno_rmFOS_NMR_imputed <- inner_join(pheno_rmHetMiss_imputed, pheno_rmFOS_NMR[,c("f.eid","FishOil")], by="f.eid")

write.table(pheno_rmFOS_NMR_imputed, "pheno_Sep2025.clean.fullPop.QCed.FOS.NMR.imputed.tsv",sep="\t",col.names=T,row.names=F,quote=F)

#################################################################################
# preprocessing oily fish intake: stratified oily fish intake into two groups: low intake group: never, less than once a week. High intake group: once a week, 2–4 times a week, 5–6 times a week, and once or more daily
#################################################################################
table(pheno_rmFOS$f.1329.0.0)
#    -3     -1      0      1      2      3      4      5
#   166   2303  48939 146689 167097  74910   3041   1025
# stratified oily fish intake into two groups: low vs high intake.
pheno_rmFOS$OilyFish <- NA
pheno_rmFOS$OilyFish[pheno_rmFOS$f.1329.0.0 >= 2] <- 1
pheno_rmFOS$OilyFish[pheno_rmFOS$f.1329.0.0 == 0 | pheno_rmFOS$f.1329.0.0 ==1] <-0
write.table(pheno_rmFOS, "pheno_Sep2025.clean.fullPop.QCed.FOS.Oily.tsv",sep="\t",col.names=T,row.names=F,quote=F)

# keep participants with NMR measurements, use omega-3 as major outcome
pheno_rmFOS_NMR <- pheno_rmFOS[!is.na(pheno_rmFOS$f.23444.0.0),]

# imputed data
pheno_rmFOS_NMR_imputed <- inner_join(pheno_rmHetMiss_imputed, pheno_rmFOS_NMR[,c("f.eid","FishOil","OilyFish")], by="f.eid")

write.table(pheno_rmFOS_NMR_imputed, "pheno_Sep2025.clean.fullPop.QCed.FOS.Oily.NMR.imputed.tsv",sep="\t",col.names=T,row.names=F,quote=F)

#################################################################################
########################### Subset each population ###############################
#################################################################################
# EUR
pheno_EUR <- pheno_rmHetMiss_imputed[pheno_rmHetMiss_imputed$pop == "EUR",c("f.eid", "f.eid")]
names(pheno_EUR) <- c("FID","IID")
write.table(pheno_EUR,"sampleID_EUR.txt",sep="\t",col.names=F,row.names=F,quote=F)

# AFR
pheno_AFR <- pheno_rmHetMiss_imputed[pheno_rmHetMiss_imputed$pop == "AFR",c("f.eid", "f.eid")]
names(pheno_AFR) <- c("FID","IID")
write.table(pheno_AFR,"sampleID_AFR.txt",sep="\t",col.names=F,row.names=F,quote=F)

# CSA
pheno_CSA <- pheno_rmHetMiss_imputed[pheno_rmHetMiss_imputed$pop == "CSA",c("f.eid", "f.eid")]
names(pheno_CSA) <- c("FID","IID")
write.table(pheno_CSA,"sampleID_CSA.txt",sep="\t",col.names=F,row.names=F,quote=F)

# EAS
pheno_EAS <- pheno_rmHetMiss_imputed[pheno_rmHetMiss_imputed$pop == "EAS",c("f.eid", "f.eid")]
names(pheno_EAS) <- c("FID","IID")
write.table(pheno_EAS,"sampleID_EAS.txt",sep="\t",col.names=F,row.names=F,quote=F)

