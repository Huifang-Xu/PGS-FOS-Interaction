'%ni%' <- Negate('%in%')
library(plyr)
library(dplyr)
library(tidyverse)
library(mice)
library(data.table)

workDir <- "phenotype/pgs"
outDir <- "phenotype/pgs"
# read in dietary information
dt <- fread(paste0(workDir, "/pheno_Sep2025.clean.fullPop.QCed.tsv"))

# list of covariates used to impute missing variables: age, sex, age × sex, body mass index, Town- send deprivation index, smoking status, alcohol status, physical activity, statin use, pop
covar <- c("f.eid", "f.31.0.0", "f.21003.0.0", "f.21001.0.0", "f.189.0.0", "f.20116.0.0", "f.20117.0.0", "f.22032.0.0", "Statin", "pop")

dt_covar <- subset(dt, select = covar)

# recode covariates: f.20116.0.0 (Smoking status), f.20117.0.0 (Alcohol drinker status), pop 
dt_covar$f.20116.0.0[dt_covar$f.20116.0.0 == -3] <- NA
dt_covar$f.20117.0.0[dt_covar$f.20117.0.0 == -3] <- NA
dt_covar$pop[dt_covar$pop == "EUR"] <- 1
dt_covar$pop[dt_covar$pop == "AFR"] <- 2
dt_covar$pop[dt_covar$pop == "AMR"] <- 3
dt_covar$pop[dt_covar$pop == "CSA"] <- 4
dt_covar$pop[dt_covar$pop == "EAS"] <- 5
dt_covar$pop[dt_covar$pop == "MID"] <- 6
dt_covar$pop <- as.integer(dt_covar$pop)

# check data structure
#'''
#for (i in 2:ncol(dt_covar)) {
#        covar_name <- names(dt_covar)[i]
#        tmp <- dt_covar %>% select(any_of(covar_name))
#        tmp <- as.data.frame(tmp)
#        print(paste0("trait: ",covar_name, "; # of unique values: ", length(unique(tmp[order(tmp[,1],decreasing=F),])))) 
#	print(unique(tmp[order(tmp[,1],decreasing=F),]))
#}
#'''

# check correlation between covariates
corr_result <- data.frame()
for (i in 2:(ncol(dt_covar)-1)) {
	for (j in (i+1):ncol(dt_covar)) {
		x_name <- names(dt_covar)[i]
		y_name <- names(dt_covar)[j]
		tmp <- subset(dt_covar, select = c(x_name,y_name))
		tmp <- tmp[complete.cases(tmp),]
		corr_test <- cor.test(tmp[[x_name]],tmp[[y_name]],method="pearson")
		corr_coef <- corr_test$estimate
		corr_pvalue <- corr_test$p.value
		corr_result <- rbind(corr_result,matrix(c(x_name,y_name,corr_coef,corr_pvalue),ncol=4))
}
}
colnames(corr_result) <- c("trait1","trait2","corr","P")
write.table(corr_result,"covar_correlation.txt",quote=F,row.names=F,col.names=T,sep="\t")

# check missing rate
total_n <- nrow(dt_covar)
for (i in 2:ncol(dt_covar)) {
	covar_name <- names(dt_covar)[i]
	missing_n <- sum(is.na(dt_covar[[covar_name]]))
	missing_rate <- missing_n/total_n
	print(paste0("trait: ",covar_name, "; missing rate: ", missing_rate))
}


# change data type and imputation methods
dt_covar <- as.data.frame(dt_covar)
for (i in 2:ncol(dt_covar)) {
	dt_covar[,i] <- as.factor(dt_covar[,i])
}
# data fields with continuous variables
continuous_variables <- c("f.21003.0.0", "f.21001.0.0", "f.189.0.0")
for (i in continuous_variables) {
	index <- which(colnames(dt_covar)==i)
	dt_covar[,index] <- as.numeric(dt_covar[,index])
}

## Impute initialization
init <- mice(dt_covar,m=5,maxit=0,seed=12345)
init$loggedEvents

## check imputation methods
## Continuous:pmm; Binary:logreg; Unordered multiple category: polyreg
impute_methods <- init$method
# replace imputation method to rf (random forest)
impute_methods <- gsub("polyreg","rf",impute_methods)
impute_methods <- gsub("pmm","rf",impute_methods)
impute_methods

## variables participate in imputation
predM <- init$predictorMatrix
predM[,1] <- 0

imputed_data <-  mice(dt_covar,m=5,maxit=20,predictorMatrix = predM, method = impute_methods,seed=12345)

saveRDS(imputed_data,file = paste0(outDir,"/pheno_Sep2025.clean.fullPop.QCed.IMPUTE_m5_maxit20.rds"))

# imputed_data <- readRDS("IMPUTE_m5_maxit20_dementia_pheno_PCs_357631.rds")

complete_data <- complete(imputed_data,1)   # Chose first imputations for further analysis

pdf(file = paste0(outDir,"/pheno_Sep2025.clean.fullPop.QCed.IMPUTE_MEAN_SD_m5_maxit20.pdf"),width=10,height = 10)
plot(imputed_data)
dev.off()
