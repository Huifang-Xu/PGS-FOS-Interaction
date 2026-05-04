'%ni%' <- Negate('%in%')
library(plyr)
library(dplyr)
library(tidyverse)
library(readr)
library(data.table)
library(rsq)
library(magrittr)
library(ggpmisc)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(RColorBrewer)
library(moments)

# accept parameters
Args <- commandArgs()

pop = Args[6]
pgs_suffix = Args[7]
phenoFile = Args[8]
inDir = Args[9]
outDir = Args[10]

# read in phenotype file, including 249 metabolites and covariates
pheno <- fread(phenoFile,header=T, sep="\t")
pheno <- as_tibble(pheno)
dim(pheno)

traits <- c("DHA", "DHA_pct", "Omega_3", "Omega_3_pct")
trait_ids <- c("f.23450.0.0", "f.23457.0.0","f.23444.0.0", "f.23451.0.0")

# output header
header_corr = c("variable1", "variable2", "pearson_corr", "pearson_pvalue", "spearman_corr","spearman_pvalue")
corr_file_to_write <- paste0(outDir,"/covariates_correlation.txt")
write.table(t(as.data.frame(header_corr)),file=corr_file_to_write,col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
m1_file_to_write <- paste0(outDir,"/asso_result_m0_24hDiet.txt")
m2_file_to_write <- paste0(outDir,"/asso_result_m2_24hDiet.txt")
m2_1_file_to_write <- paste0(outDir,"/asso_result_m2_1_24hDiet.txt")
m2_3_file_to_write <- paste0(outDir,"/asso_result_m2_3_24hDiet.txt")


#####################################################################
########################## Various functions ########################
#####################################################################
# rank-based inverse transformation
inversenormal <- function(x) {
    # inverse normal if you have missing data
    return(qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))
}

# function of removing outlier in data using interquantile range
outlierRemoval <- function(data, col_name, k=1.5) {
        qs  <- quantile(data[[col_name]], c(.25, .75), na.rm = TRUE, type = 7)
        iqr <- diff(qs)
        low <- qs[1] - k*iqr
        high <- qs[2] + k*iqr
        out_data <- data[data[[col_name]] >= low & data[[col_name]] <= high,]
        return(out_data)
}

# function of scatter plots for FOS
colors <- brewer.pal(5,"Set1")[1:2]
plot_scatter1 <- function(png_name, data,title,x_name,y_name, width, height, beta, pvalue) {
#  pdf(pdf_name, width = pdf_width, height=pdf_height)
png(png_name, type=c("cairo", "cairo-png", "Xlib", "quartz"), width=width, height=height, units="in", res=720)
  plt <- ggplot(data=data,aes(x=x,y=y,colour=group))+
    geom_point(data=data,size=1.3,aes(x=x,y=y,colour=group), alpha=0.3)+
    stat_poly_line(data=data[data$group=="FOS",],aes(x=x,y=y),colour=colors[1],linewidth=0.8,formula = y~x) +
    stat_poly_eq(data = data[data$group=="FOS",],aes(x=x,y=y,label = paste(after_stat(eq.label),after_stat(rr.label),"(FOS)", sep = "~~~")),formula = y~x, parse = TRUE,label.y = 0.98,colour=colors[1])+ #y~x+0
    # two groups
	stat_poly_line(data=data[data$group=="non-FOS",],aes(x=x,y=y),colour=colors[2],linewidth=0.8,formula = y~x)+
	stat_poly_eq(data = data[data$group=="non-FOS",],aes(x=x,y=y,label = paste(after_stat(eq.label),after_stat(rr.label),"(non-FOS)",sep = "~~~")),formula = y~x, parse = TRUE,label.y = 0.93,colour=colors[2])+
    theme_bw()+
    theme(legend.title=element_blank(),
          legend.position=c(0.1,0.1),
          legend.background = element_rect(fill = 'white', colour = 'black'),
          axis.title.x=element_text(size=13),
          axis.title.y=element_text(size=13),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12))+
    labs(title = title)+
    xlab(x_name) +
    ylab(y_name) +
    scale_color_manual(values=colors)
  print(plt)
  dev.off()
}

plot_scatter2 <- function(png_name, data,title,x_name,y_name, width, height, beta, pvalue) {
#  pdf(pdf_name, width = pdf_width, height=pdf_height)
png(png_name, type=c("cairo", "cairo-png", "Xlib", "quartz"), width=width, height=height, units="in", res=720)
  plt <- ggplot(data=data,aes(x=x,y=y,colour=group))+
    geom_point(data=data,size=1.3,aes(x=x,y=y,colour=group), alpha=0.3)+
    stat_poly_line(data=data[data$group=="high-OFI",],aes(x=x,y=y),colour=colors[1],linewidth=0.8,formula = y~x) +
    stat_poly_eq(data = data[data$group=="high-OFI",],aes(x=x,y=y,label = paste(after_stat(eq.label),after_stat(rr.label),"(high-OFI)", sep = "~~~")),formula = y~x, parse = TRUE,label.y = 0.98,colour=colors[1])+ #y~x+0
    # two groups
        stat_poly_line(data=data[data$group=="low-OFI",],aes(x=x,y=y),colour=colors[2],linewidth=0.8,formula = y~x)+
        stat_poly_eq(data = data[data$group=="low-OFI",],aes(x=x,y=y,label = paste(after_stat(eq.label),after_stat(rr.label),"(low-OFI)",sep = "~~~")),formula = y~x, parse = TRUE,label.y = 0.93,colour=colors[2])+
    theme_bw()+
    theme(legend.title=element_blank(),
          legend.position=c(0.1,0.1),
          legend.background = element_rect(fill = 'white', colour = 'black'),
          axis.title.x=element_text(size=13),
          axis.title.y=element_text(size=13),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12))+
    labs(title = title)+
    xlab(x_name) +
    ylab(y_name) +
    scale_color_manual(values=colors)
  print(plt)
  dev.off()
}

# function of boxplot
colors <- brewer.pal(5,"Set1")[1:2]
plot_boxplot <- function(png_name,data,title,y_name, width, height) {
png(png_name, type=c("cairo", "cairo-png", "Xlib", "quartz"), width=width, height=height, units="in", res=720)
#  pdf(pdf_name, width = pdf_width, height=pdf_height)
#  plt <- ggplot(data=data,aes(x=group,y=y,fill=FOS_24h_intake,group = interaction(group,FOS_24h_intake)))+
plt <- ggplot(data=data,aes(x=pgs_group,y=y,fill=fo_group,group = interaction(pgs_group,fo_group)))+
	geom_violin(alpha=0.4,position=position_dodge(0.70))+
	  theme_bw()+
    theme(legend.title=element_blank(),
          legend.position="top",
          legend.background = element_rect(fill = 'white', colour = 'black'),
          axis.title.x=element_text(size=13),
          axis.title.y=element_text(size=13),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12))+
    labs(title = title)+
    xlab("PGS percentile") +
    ylab(y_name) +
    geom_boxplot(width=0.1,fill="white",outlier.colour = NA,position=position_dodge(0.70))+ 
    stat_summary(fun = mean,geom="point",fill="black",shape=21,size=2.5,position=position_dodge(0.70))
    scale_fill_manual(values=colors)
  print(plt)
  dev.off()
}

# function of testing correlation between two variables
corr_test <- function(data_in, x_name, y_name, file_to_write) {
	data_xy <- subset(data_in, select = c(x_name, y_name))
	data_xy <- data_xy[complete.cases(data_xy),]
	pearson_corr <- cor.test(x= data_xy[[x_name]], y = data_xy[[y_name]], method = "pearson")$estimate
        pearson_pvalue <- format.pval(cor.test(x= data_xy[[x_name]], y = data_xy[[y_name]], method = "pearson")$p.value, digits = 10)
        spearman_corr <- cor.test(x= data_xy[[x_name]], y = data_xy[[y_name]], method = "spearman")$estimate
        spearman_pvalue <- format.pval(cor.test(x= data_xy[[x_name]], y = data_xy[[y_name]], method = "spearman")$p.value, digits = 10)
	corr_result <- t(as.data.frame(c(x_name, y_name, pearson_corr, pearson_pvalue, spearman_corr, spearman_pvalue)))
        write.table(corr_result, file=file_to_write,col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
}


assc_test <- function(data_assc,  y_name, formula, group, file_to_write) {
	data_assc <- data_assc[complete.cases(data_assc),]
	fitModel <- lm(formula, data=data_assc)
	sampleSize <- nobs(fitModel)
        ml_summary <- summary(fitModel)
        model_rsq <- rsq(fitModel)
        rsqValue <- rsq.partial(fitModel, adj = TRUE)
        # convert lm results into one row
	coefs <- as.data.frame(ml_summary$coefficients, check.names = FALSE)
	colnames(coefs) <- c("coef", "se", "tStats", "pvalue")
	rownames(coefs)[rownames(coefs)=='(Intercept)'] <- 'Intercept'
	coef_results <- coefs %>%  tibble::rownames_to_column("term") %>%  pivot_longer(-term, names_to = "stat", values_to = "value") %>%  mutate(name = paste(term, stat, sep = "_")) %>%  select(name, value) %>%  pivot_wider(names_from = name, values_from = value)
        # convert rsq results into one row
	rsq_vars <- paste0(rsqValue$variable, "_rsq")
  	rsq_values <- rsqValue$partial.rsq
	rsq_results <- as.data.frame(setNames(as.list(rsq_values), rsq_vars), check.names = FALSE)
	# output final result
	result <- bind_cols(tibble(sampleSize = sampleSize, trait_name  = trait_name, group = group), coef_results,tibble(model_rsq = model_rsq),rsq_results)
	write.table(result, file=file_to_write,col.names = TRUE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
}

#####################################################################
#####################################################################
#####################################################################


#####################################################################
############### data processing and association test ################
#####################################################################
# subset target population
pheno = pheno[pheno$pop == pop,]
# subset covariates
covars_all <- c("f.31.0.0", "f.21003.0.0", "sex_by_age", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "BMI_imp", "TSI_imp", "smoking_imp", "alcohol_imp", "physical_imp", "Statin", "FOS_24h")
pheno_covar <- subset(pheno, select = c("f.eid", trait_ids,covars_all))

# test correlation between covariates
for (c in 1:(length(covars_all)-1)) {
        for (i in (c+1):length(covars_all)) {
                # test correlation between covariates
 corr_test(pheno_covar, covars_all[c], covars_all[i], corr_file_to_write)
        }
}

# iterate each trait
for (t in 1:length(traits)) {
trait_name = traits[t]
trait_field = trait_ids[t]

# read in pgs file
pgs <- fread(paste0(inDir,"/",trait_name,pgs_suffix))


# merge raw phenotype with pgs
pheno.pgs <- inner_join(pheno_covar,pgs[,c("IID","SCORE")],by=c("f.eid"="IID"))

# remove outliers in raw data
pheno.pgs.rmOutli <- outlierRemoval(pheno.pgs, trait_field, k=1.5)

# scaled PGS
pheno.pgs.rmOutli$SCORE.norm <- as.vector(scale(pheno.pgs.rmOutli$SCORE, center = TRUE, scale = TRUE))
# add FOS-by-PGS interaction term
pheno.pgs.rmOutli$FOS_PGS <- pheno.pgs.rmOutli$FOS_24h * pheno.pgs.rmOutli$SCORE.norm

pheno.pgs.raw <- pheno.pgs
pheno.pgs <- pheno.pgs.rmOutli

# rank-based inverse normalization
pheno.pgs.INT <- copy(pheno.pgs)
pheno.pgs.INT[[trait_field]] <- inversenormal(pheno.pgs.INT[[trait_field]])

# test correlation between interaction term and covariats
for (c in 1:length(covars_all)) {
	corr_test(pheno.pgs, "FOS_PGS", covars_all[c], corr_file_to_write)
	corr_test(pheno.pgs, "OFI_PGS", covars_all[c], corr_file_to_write)
}
# test correlation between interaction term and phenotype
dt_to_corr <- subset(pheno.pgs, select = c(trait_field,"FOS_PGS","OFI_PGS"))
names(dt_to_corr) <- c(paste0(trait_field,"_raw"),"FOS_PGS","OFI_PGS")
corr_test(dt_to_corr, paste0(trait_field,"_raw"), "FOS_PGS",corr_file_to_write)
corr_test(dt_to_corr, paste0(trait_field,"_raw"), "OFI_PGS",corr_file_to_write)
dt_to_corr <- subset(pheno.pgs.INT, select = c(trait_field,"FOS_PGS", "OFI_PGS"))
names(dt_to_corr) <- c(paste0(trait_field,"_INT"),"FOS_PGS", "OFI_PGS")
corr_test(dt_to_corr, paste0(trait_field,"_INT"), "FOS_PGS",corr_file_to_write)
corr_test(dt_to_corr, paste0(trait_field,"_INT"), "OFI_PGS",corr_file_to_write)
# test correlation between interaction term and PGS
corr_test(pheno.pgs, "SCORE.norm","FOS_PGS", corr_file_to_write)
corr_test(pheno.pgs, "SCORE.norm","OFI_PGS", corr_file_to_write)

# test correlation between raw and INT phenotype
dt_to_corr <- subset(pheno.pgs, select = c("f.eid", trait_field))
names(dt_to_corr) <- c("f.eid", paste0(trait_field,"_raw"))
dt_to_corr <- left_join(dt_to_corr, subset(pheno.pgs.INT,select = c("f.eid", trait_field)), by = "f.eid")
colnames(dt_to_corr) <- c("f.eid", paste0(trait_field,"_raw"),paste0(trait_field,"_INT"))
corr_test(dt_to_corr, paste0(trait_field,"_raw"), paste0(trait_field,"_INT"), corr_file_to_write)

# test correlation between PGS and phenotypes
dt_to_corr <- subset(pheno.pgs, select = c("f.eid", trait_field,"SCORE.norm"))
names(dt_to_corr) <- c("f.eid", paste0(trait_field,"_raw"),paste0("PGS_",trait_field))
corr_test(dt_to_corr,paste0(trait_field,"_raw"),paste0("PGS_",trait_field), corr_file_to_write)
dt_to_corr <- subset(pheno.pgs.INT, select = c("f.eid", trait_field,"SCORE.norm"))
names(dt_to_corr) <- c("f.eid", paste0(trait_field,"_INT"),paste0("PGS_",trait_field))
corr_test(dt_to_corr,paste0(trait_field,"_INT"),paste0("PGS_",trait_field), corr_file_to_write)


# test correlation between phenotype/PGS and covariates
for (c in 1:length(covars_all)) {
	# test correlation between raw phenotype and covariates
	dt_to_corr <- subset(pheno.pgs, select = c("f.eid", trait_field,covars_all[c]))
	names(dt_to_corr) <- c("f.eid", paste0(trait_field,"_raw"),covars_all[c])
	corr_test(dt_to_corr,paste0(trait_field,"_raw"),covars_all[c],corr_file_to_write)
	# test correlation between transformed phenotype and covariates
	dt_to_corr <- subset(pheno.pgs.INT, select = c("f.eid", trait_field,covars_all[c]))
	names(dt_to_corr) <- c("f.eid", paste0(trait_field,"_INT"),covars_all[c])
        corr_test(dt_to_corr,paste0(trait_field,"_INT"),covars_all[c],corr_file_to_write)
	# test correlation between each PGS and covariates
	dt_to_corr <- subset(pheno.pgs,select = c("f.eid","SCORE.norm",covars_all[c]))
	colnames(dt_to_corr) <- c("f.eid",paste0("PGS_",trait_field), covars_all[c])
	corr_test(dt_to_corr,paste0("PGS_",trait_field),covars_all[c],corr_file_to_write)
}


###########################
# strategy 1.0: stratify by two groups: fish oil consumer yes vs no; all participants
dt_s1 <- copy(pheno.pgs)
dt_s1_g1 <- dt_s1[dt_s1$FOS_24h==1,]
dt_s1_g2 <- dt_s1[dt_s1$FOS_24h==0,]
# rank-based inverse normalization
dt_s1_INT <- copy(pheno.pgs.INT)
dt_s1_g1_INT <- dt_s1_INT[dt_s1_INT$FOS_24h==1,]
dt_s1_g2_INT <- dt_s1_INT[dt_s1_INT$FOS_24h==0,]

###########################

###########################
# strategy 2: stratify by two groups: participants with top 5% vs bottom 5% PGS
dt_s2 <- copy(pheno.pgs)
# Create percentile breaks
dt_s2$pgs_perct <- cut(
  rank(dt_s2$SCORE.norm) / length(dt_s2$SCORE.norm),
  breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25,0.3, 0.35, 0.4,0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1),
  labels = c("0-5%", "5-10%","10-15%","15-20%","20-25%","25-30%","30-35%","35-40%", "40-45%","45-50%","50-55%","55-60%","60-65%", "65-70%","70-75%","75-80%","80-85%","85-90%","90-95%","95-100%"),
  include.lowest = TRUE
)
# subset participants in 0-5% PGS percentile
dt_s2_g1 <- dt_s2[dt_s2$pgs_perct == "0-5%",]
# subset participants in 95-100% PGS percentile
dt_s2_g2 <- dt_s2[dt_s2$pgs_perct == "95-100%",]

# rank-based inverse normalization
dt_s2_INT <- copy(pheno.pgs.INT)
dt_s2_INT$pgs_perct <- cut(
  rank(dt_s2_INT$SCORE.norm) / length(dt_s2_INT$SCORE.norm),
  breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25,0.3, 0.35, 0.4,0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1),
  labels = c("0-5%", "5-10%","10-15%","15-20%","20-25%","25-30%","30-35%","35-40%", "40-45%","45-50%","50-55%","55-60%","60-65%", "65-70%","70-75%","75-80%","80-85%","85-90%","90-95%","95-100%"),
  include.lowest = TRUE
)
dt_s2_g1_INT <- dt_s2_INT[dt_s2_INT$pgs_perct == "0-5%",]
dt_s2_g2_INT <- dt_s2_INT[dt_s2_INT$pgs_perct == "95-100%",]

###########################

###########################
# strategy 3: stratify by five groups: participants in 0-20%, 21-40%, 41-60%, 61-80%, 81-100% PGS percentile
dt_s3 <- copy(pheno.pgs)
# Create percentile breaks
dt_s3$pgs_perct <- cut(
  rank(dt_s3$SCORE.norm) / length(dt_s3$SCORE.norm),
  breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
  labels = c("0-20%", "20-40%", "40-60%", "60-80%", "80-100%"),
  include.lowest = TRUE
)
# subset 5 groups
dt_s3_g1 <- dt_s3[dt_s3$pgs_perct =="0-20%",]
dt_s3_g2 <- dt_s3[dt_s3$pgs_perct =="20-40%",]
dt_s3_g3 <- dt_s3[dt_s3$pgs_perct =="40-60%",]
dt_s3_g4 <- dt_s3[dt_s3$pgs_perct =="60-80%",]
dt_s3_g5 <- dt_s3[dt_s3$pgs_perct =="80-100%",]

# rank-based inverse normalization
dt_s3_INT <- copy(pheno.pgs.INT)
# Create percentile breaks
dt_s3_INT$pgs_perct <- cut(
  rank(dt_s3_INT$SCORE.norm) / length(dt_s3_INT$SCORE.norm),
  breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
  labels = c("0-20%", "20-40%", "40-60%", "60-80%", "80-100%"),
  include.lowest = TRUE
)
# subset 5 groups
dt_s3_g1_INT <- dt_s3_INT[dt_s3_INT$pgs_perct =="0-20%",]
dt_s3_g2_INT <- dt_s3_INT[dt_s3_INT$pgs_perct =="20-40%",]
dt_s3_g3_INT <- dt_s3_INT[dt_s3_INT$pgs_perct =="40-60%",]
dt_s3_g4_INT <- dt_s3_INT[dt_s3_INT$pgs_perct =="60-80%",]
dt_s3_g5_INT <- dt_s3_INT[dt_s3_INT$pgs_perct =="80-100%",]

###########################

###########################
# construct formula of association test
# model1 (w/o interaction): PUFA ~ sex (f.31.0.0) + age (f.21003.0.0) + sex * age + BMI + TSI + smoking + alcohol + physical actvity + statin use + FOS
covar_m1 = c("f.31.0.0","f.21003.0.0","sex_by_age", "BMI_imp", "TSI_imp", "smoking_imp", "alcohol_imp", "physical_imp", "Statin", "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","FOS_24h", "SCORE.norm")
formula_m1 <- as.formula(paste(trait_field, "~", paste(covar_m1, collapse = " + ")))

# model2 (w/ interaction): PUFA ~ sex (f.31.0.0) + age (f.21003.0.0) + PC1 ~ PC10 + sex * age + BMI + TSI + smoking + alcohol + physical actvity + statin use + FOS + PGS (norm) + FOS * PGS
covar_m2 = c("f.31.0.0","f.21003.0.0","sex_by_age","BMI_imp", "TSI_imp", "smoking_imp", "alcohol_imp", "physical_imp", "Statin", "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10", "FOS_24h","SCORE.norm","FOS_PGS")
formula_m2 <- as.formula(paste(trait_field, "~", paste(covar_m2, collapse = " + ")))

# model2.1 (stratified by FOS): PUFA ~ sex (f.31.0.0) + age (f.21003.0.0) + PC1 ~ PC10 + sex * age + BMI + TSI + smoking + alcohol + physical actvity + statin use + PGS (norm) 
covar_m2_1 = c("f.31.0.0","f.21003.0.0","sex_by_age","BMI_imp", "TSI_imp", "smoking_imp", "alcohol_imp", "physical_imp", "Statin", "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10", "SCORE.norm")
formula_m2_1 <- as.formula(paste(trait_field, "~", paste(covar_m2_1, collapse = " + ")))

# model2.3 (stratified by PGS, without consider FOS*PGS interaction term): PUFA ~ sex (f.31.0.0) + age (f.21003.0.0) + PC1 ~ PC10 + sex * age + BMI + TSI + smoking + alcohol + physical actvity + statin use + FOS + PGS (norm) 
covar_m2_3 = c("f.31.0.0","f.21003.0.0","sex_by_age","BMI_imp", "TSI_imp", "smoking_imp", "alcohol_imp", "physical_imp", "Statin", "FOS_24h")
formula_m2_3 <- as.formula(paste(trait_field, "~", paste(covar_m2_3, collapse = " + ")))

###########################

###########################
######### perform association test: full set; model1 - model2; raw and INT phenotype
# model: m1; phenotype: raw; primary analysis
data_assc <- subset(pheno.pgs, select=c(trait_field,covar_m1))
assc_test(data_assc, y_name = trait_field, formula = formula_m1, group = paste0(pop,"_all_m1_raw"), file_to_write = m1_file_to_write)
# model: m1; phenotype: rank-based inverse normalization; primary analysis
data_assc <- subset(pheno.pgs.INT, select=c(trait_field,covar_m1))
assc_test(data_assc,  y_name = trait_field, formula = formula_m1, group = paste0(pop,"_all_m1_INT"), file_to_write = m1_file_to_write)

# model: m2; phenotype: raw; primary analysis
data_assc <- subset(pheno.pgs, select=c(trait_field,covar_m2))
assc_test(data_assc,  y_name = trait_field, formula = formula_m2, group = paste0(pop,"_all_m2_raw_primary"), file_to_write = m2_file_to_write)
# model: m2; phenotype: rank-based inverse normalization; primary analysis
data_assc <- subset(pheno.pgs.INT, select=c(trait_field,covar_m2))
assc_test(data_assc, y_name = trait_field, formula = formula_m2, group = paste0(pop,"_all_m2_INT_primary"), file_to_write = m2_file_to_write)

############ perform association test: subgroup (stratified by FOS, or PGS groups); model2.1 - model 2.3; raw and INT phenotype
# model: m2.1; phenotype: raw; stratified by FOS_24h
data_assc <- subset(dt_s1_g1, select=c(trait_field,covar_m2_1))
assc_test(data_assc,  y_name = trait_field, formula = formula_m2_1, group = paste0(pop,"_FOS_m21_raw"), file_to_write = m2_1_file_to_write)
data_assc <- subset(dt_s1_g2, select=c(trait_field,covar_m2_1))
assc_test(data_assc,  y_name = trait_field, formula = formula_m2_1, group = paste0(pop,"_nonFOS_m21_raw"), file_to_write = m2_1_file_to_write)
# model: m21; phenotype: rank-based inverse normalization
data_assc <- subset(dt_s1_g1_INT, select=c(trait_field,covar_m2_1))
assc_test(data_assc,  y_name = trait_field, formula = formula_m2_1, group = paste0(pop,"_FOS_m21_INT"), file_to_write = m2_1_file_to_write)
data_assc <- subset(dt_s1_g2_INT, select=c(trait_field,covar_m2_1))
assc_test(data_assc,  y_name = trait_field, formula = formula_m2_1, group = paste0(pop,"_nonFOS_m21_INT"), file_to_write = m2_1_file_to_write)

# model: m2.3; phenotype: raw; all PGS 
data_assc <- subset(pheno.pgs, select=c(trait_field,covar_m2_3))
assc_test(data_assc,  y_name = trait_field, formula = formula_m2_3, group = paste0(pop,"_allPGS_m23_raw"), file_to_write = m2_3_file_to_write)
# model: m23; phenotype: rank-based inverse normalization; allPGS
data_assc <- subset(pheno.pgs.INT, select=c(trait_field,covar_m2_3))
assc_test(data_assc,  y_name = trait_field, formula = formula_m2_3, group = paste0(pop,"_allPGS_m23_INT"), file_to_write = m2_3_file_to_write)
# model: m23; phenotype: raw; stratified by PGS bottom 5% vs top 5%
data_assc <- subset(dt_s2_g1, select=c(trait_field,covar_m2_3))
assc_test(data_assc,  y_name = trait_field, formula = formula_m2_3, group = paste0(pop,"_5%PGS_m23_raw"), file_to_write = m2_3_file_to_write)
data_assc <- subset(dt_s2_g2, select=c(trait_field,covar_m2_3))
assc_test(data_assc,  y_name = trait_field, formula = formula_m2_3, group = paste0(pop,"_95%PGS_m23_raw"), file_to_write = m2_3_file_to_write)
# model: m2.3; phenotype: rank-based inverse normalization; stratified by PGS bottom 5% vs top 5%
data_assc <- subset(dt_s2_g1_INT, select=c(trait_field,covar_m2_3))
assc_test(data_assc,  y_name = trait_field, formula = formula_m2_3, group = paste0(pop,"_5%PGS_m23_INT"), file_to_write = m2_3_file_to_write)
data_assc <- subset(dt_s2_g2_INT, select=c(trait_field,covar_m2_3))
assc_test(data_assc,  y_name = trait_field, formula = formula_m2_3, group = paste0(pop,"_95%PGS_m23_INT"), file_to_write = m2_3_file_to_write)

# model: m2.3; phenotype: raw; stratified by PGS 0-20, 20-40, 40-60, 60-80, 80-100%
data_assc <- subset(dt_s3_g1, select=c(trait_field,covar_m2_3))
assc_test(data_assc,  y_name = trait_field, formula = formula_m2_3, group = paste0(pop,"_0-20%PGS_m23_raw"), file_to_write = m2_3_file_to_write)
data_assc <- subset(dt_s3_g2, select=c(trait_field,covar_m2_3))
assc_test(data_assc,  y_name = trait_field, formula = formula_m2_3, group = paste0(pop,"_20-40%PGS_m23_raw"), file_to_write = m2_3_file_to_write)
data_assc <- subset(dt_s3_g3, select=c(trait_field,covar_m2_3))
assc_test(data_assc,  y_name = trait_field, formula = formula_m2_3, group = paste0(pop,"_40-60%PGS_m23_raw"), file_to_write = m2_3_file_to_write)
data_assc <- subset(dt_s3_g4, select=c(trait_field,covar_m2_3))
assc_test(data_assc,  y_name = trait_field, formula = formula_m2_3, group = paste0(pop,"_60-80%PGS_m23_raw"), file_to_write = m2_3_file_to_write)
data_assc <- subset(dt_s3_g5, select=c(trait_field,covar_m2_3))
assc_test(data_assc,  y_name = trait_field, formula = formula_m2_3, group = paste0(pop,"_80-100%PGS_m23_raw"), file_to_write = m2_3_file_to_write)
# model: m2.3; phenotype: rank-based inverse normalization; stratified by PGS 0-20, 20-40, 40-60, 60-80, 80-100%
data_assc <- subset(dt_s3_g1_INT, select=c(trait_field,covar_m2_3))
assc_test(data_assc,  y_name = trait_field, formula = formula_m2_3, group = paste0(pop,"_0-20%PGS_m23_INT"), file_to_write = m2_3_file_to_write)
data_assc <- subset(dt_s3_g2_INT, select=c(trait_field,covar_m2_3))
assc_test(data_assc,  y_name = trait_field, formula = formula_m2_3, group = paste0(pop,"_20-40%PGS_m23_INT"), file_to_write = m2_3_file_to_write)
data_assc <- subset(dt_s3_g3_INT, select=c(trait_field,covar_m2_3))
assc_test(data_assc,  y_name = trait_field, formula = formula_m2_3, group = paste0(pop,"_40-60%PGS_m23_INT"), file_to_write = m2_3_file_to_write)
data_assc <- subset(dt_s3_g4_INT, select=c(trait_field,covar_m2_3))
assc_test(data_assc,  y_name = trait_field, formula = formula_m2_3, group = paste0(pop,"_60-80%PGS_m23_INT"), file_to_write = m2_3_file_to_write)
data_assc <- subset(dt_s3_g5_INT, select=c(trait_field,covar_m2_3))
assc_test(data_assc,  y_name = trait_field, formula = formula_m2_3, group = paste0(pop,"_80-100%PGS_m23_INT"), file_to_write = m2_3_file_to_write)

###########################
# plot

# plot scatter plot: all PGS; Fish oil
df_to_plot <- subset(dt_s1, select=c("SCORE.norm", trait_field, "FOS_24h"))
colnames(df_to_plot) <- c("x", "y", "group")
df_to_plot$group[df_to_plot$group==1] <- "FOS"
df_to_plot$group[df_to_plot$group==0] <- "non-FOS"
plot_scatter1(png_name=paste0(outDir,"/scatter_strat_allPGS_FOS_",trait_name,".png"),data = df_to_plot,title = paste0(trait_name," (raw); all PGS groups; FOS"),x_name = 'PGS (scaled)',y_name = trait_name,width = 6, height = 5)

# plot scatter plot: bottom 5% PGS; Fish oil
df_to_plot <- subset(dt_s2_g1, select=c("SCORE.norm", trait_field, "FOS_24h"))
colnames(df_to_plot) <- c("x", "y", "group")
df_to_plot$group[df_to_plot$group==1] <- "FOS"
df_to_plot$group[df_to_plot$group==0] <- "non-FOS"
plot_scatter1(png_name=paste0(outDir,"/scatter_strat_bottom5PGS_FOS_",trait_name,".png"),data = df_to_plot,title = paste0(trait_name," (raw); bottom 5% PGS groups; FOS"),x_name = 'PGS (scaled)',y_name = trait_name,width = 6, height = 5)

# plot scatter plot: top 5% PGS; Fish oil
df_to_plot <- subset(dt_s2_g2, select=c("SCORE.norm", trait_field, "FOS_24h"))
colnames(df_to_plot) <- c("x", "y", "group")
df_to_plot$group[df_to_plot$group==1] <- "FOS"
df_to_plot$group[df_to_plot$group==0] <- "non-FOS"
plot_scatter1(png_name=paste0(outDir,"/scatter_strat_top5PGS_FOS_",trait_name,".png"),data = df_to_plot,title = paste0(trait_name," (raw); top 5% PGS groups; FOS"),x_name = 'PGS (scaled)',y_name = trait_name,width = 6, height = 5)

# plot scatter plot: 0-20% PGS; Fish oil
df_to_plot <- subset(dt_s3_g1, select=c("SCORE.norm", trait_field, "FOS_24h"))
colnames(df_to_plot) <- c("x", "y", "group")
df_to_plot$group[df_to_plot$group==1] <- "FOS"
df_to_plot$group[df_to_plot$group==0] <- "non-FOS"
plot_scatter1(png_name=paste0(outDir,"/scatter_strat_0To20PGS_FOS_",trait_name,".png"),data = df_to_plot,title = paste0(trait_name," (raw); 0-20% PGS groups; FOS"),x_name = 'PGS (scaled)',y_name = trait_name,width = 6, height = 5)

# plot scatter plot: 21-40% PGS: Fish oil
df_to_plot <- subset(dt_s3_g2, select=c("SCORE.norm", trait_field, "FOS_24h"))
colnames(df_to_plot) <- c("x", "y", "group")
df_to_plot$group[df_to_plot$group==1] <- "FOS"
df_to_plot$group[df_to_plot$group==0] <- "non-FOS"
plot_scatter1(png_name=paste0(outDir,"/scatter_strat_20To40PGS_FOS_",trait_name,".png"),data = df_to_plot,title = paste0(trait_name," (raw); 20-40% PGS groups; FOS"),x_name = 'PGS (scaled)',y_name = trait_name,width = 6, height = 5)

# plot scatter plot: 41-60% PGS; Fish oil
df_to_plot <- subset(dt_s3_g3, select=c("SCORE.norm", trait_field, "FOS_24h"))
colnames(df_to_plot) <- c("x", "y", "group")
df_to_plot$group[df_to_plot$group==1] <- "FOS"
df_to_plot$group[df_to_plot$group==0] <- "non-FOS"
plot_scatter1(png_name=paste0(outDir,"/scatter_strat_40To60PGS_FOS_",trait_name,".png"),data = df_to_plot,title = paste0(trait_name," (raw); 40-60% PGS groups; FOS"),x_name = 'PGS (scaled)',y_name = trait_name,width = 6, height = 5)

# plot scatter plot: 61-80% PGS; Fish oil
df_to_plot <- subset(dt_s3_g4, select=c("SCORE.norm", trait_field, "FOS_24h"))
colnames(df_to_plot) <- c("x", "y", "group")
df_to_plot$group[df_to_plot$group==1] <- "FOS"
df_to_plot$group[df_to_plot$group==0] <- "non-FOS"
plot_scatter1(png_name=paste0(outDir,"/scatter_strat_60To80PGS_FOS_",trait_name,".png"),data = df_to_plot,title = paste0(trait_name," (raw); 60-80% PGS groups; FOS"),x_name = 'PGS (scaled)',y_name = trait_name,width = 6, height = 5)

# plot scatter plot: 81-100% PGS; Fish oil
df_to_plot <- subset(dt_s3_g5, select=c("SCORE.norm", trait_field, "FOS_24h"))
colnames(df_to_plot) <- c("x", "y", "group")
df_to_plot$group[df_to_plot$group==1] <- "FOS"
df_to_plot$group[df_to_plot$group==0] <- "non-FOS"
plot_scatter1(png_name=paste0(outDir,"/scatter_strat_80To100PGS_FOS_",trait_name,".png"),data = df_to_plot,title = paste0(trait_name," (raw); 80-100% PGS groups; FOS"),x_name = 'PGS (scaled)',y_name = trait_name,width = 6, height = 5)


# plot scatter plot: all PGS; INT normalized phenotype; FOS_24h
df_to_plot <- subset(dt_s1_INT, select=c("SCORE.norm", trait_field, "FOS_24h"))
colnames(df_to_plot) <- c("x", "y", "group")
df_to_plot$group[df_to_plot$group==1] <- "FOS"
df_to_plot$group[df_to_plot$group==0] <- "non-FOS"
plot_scatter1(png_name=paste0(outDir,"/scatter_strat_allPGS_INT_FOS_",trait_name,".png"),data = df_to_plot,title = paste0(trait_name," (INT); all PGS groups; FOS"),x_name = 'PGS (scaled)',y_name = trait_name,width = 6, height = 5)

# plot scatter plot: bottom 5% PGS; INT normalized phenotype; FOS_24h
df_to_plot <- subset(dt_s2_g1_INT, select=c("SCORE.norm", trait_field, "FOS_24h"))
colnames(df_to_plot) <- c("x", "y", "group")
df_to_plot$group[df_to_plot$group==1] <- "FOS"
df_to_plot$group[df_to_plot$group==0] <- "non-FOS"
plot_scatter1(png_name=paste0(outDir,"/scatter_strat_bottom5PGS_INT_FOS_",trait_name,".png"),data = df_to_plot,title = paste0(trait_name," (INT); bottom 5% PGS groups; FOS"),x_name = 'PGS (scaled)',y_name = trait_name,width = 6, height = 5)

# plot scatter plot: top 5% PGS; INT normalized phenotype; FOS_24h
df_to_plot <- subset(dt_s2_g2_INT, select=c("SCORE.norm", trait_field, "FOS_24h"))
colnames(df_to_plot) <- c("x", "y", "group")
df_to_plot$group[df_to_plot$group==1] <- "FOS"
df_to_plot$group[df_to_plot$group==0] <- "non-FOS"
plot_scatter1(png_name=paste0(outDir,"/scatter_strat_top5PGS_INT_FOS_",trait_name,".png"),data = df_to_plot,title = paste0(trait_name," (INT); top 5% PGS groups; FOS"),x_name = 'PGS (scaled)',y_name = trait_name,width = 6, height = 5)

# plot scatter plot: 0-20% PGS; INT normalized phenotype; FOS_24h
df_to_plot <- subset(dt_s3_g1_INT, select=c("SCORE.norm", trait_field, "FOS_24h"))
colnames(df_to_plot) <- c("x", "y", "group")
df_to_plot$group[df_to_plot$group==1] <- "FOS"
df_to_plot$group[df_to_plot$group==0] <- "non-FOS"
plot_scatter1(png_name=paste0(outDir,"/scatter_strat_0To20PGS_INT_FOS_",trait_name,".png"),data = df_to_plot,title = paste0(trait_name," (INT); 0-20% PGS groups; FOS"),x_name = 'PGS (scaled)',y_name = trait_name,width = 6, height = 5)

# plot scatter plot: 21-40% PGS; INT normalized phenotype; FOS_24h
df_to_plot <- subset(dt_s3_g2_INT, select=c("SCORE.norm", trait_field, "FOS_24h"))
colnames(df_to_plot) <- c("x", "y", "group")
df_to_plot$group[df_to_plot$group==1] <- "FOS"
df_to_plot$group[df_to_plot$group==0] <- "non-FOS"
plot_scatter1(png_name=paste0(outDir,"/scatter_strat_20To40PGS_INT_FOS_",trait_name,".png"),data = df_to_plot,title = paste0(trait_name," (INT); 20-40% PGS groups; FOS"),x_name = 'PGS (scaled)',y_name = trait_name,width = 6, height = 5)

# plot scatter plot: 41-60% PGS; INT normalized phenotype; FOS_24h
df_to_plot <- subset(dt_s3_g3_INT, select=c("SCORE.norm", trait_field, "FOS_24h"))
colnames(df_to_plot) <- c("x", "y", "group")
df_to_plot$group[df_to_plot$group==1] <- "FOS"
df_to_plot$group[df_to_plot$group==0] <- "non-FOS"
plot_scatter1(png_name=paste0(outDir,"/scatter_strat_40To60PGS_INT_FOS_",trait_name,".png"),data = df_to_plot,title = paste0(trait_name," (INT); 40-60% PGS groups; FOS"),x_name = 'PGS (scaled)',y_name = trait_name,width = 6, height = 5)

# plot scatter plot: 61-80% PGS; INT normalized phenotype; FOS_24h
df_to_plot <- subset(dt_s3_g4_INT, select=c("SCORE.norm", trait_field, "FOS_24h"))
colnames(df_to_plot) <- c("x", "y", "group")
df_to_plot$group[df_to_plot$group==1] <- "FOS"
df_to_plot$group[df_to_plot$group==0] <- "non-FOS"
plot_scatter1(png_name=paste0(outDir,"/scatter_strat_60To80PGS_INT_FOS_",trait_name,".png"),data = df_to_plot,title = paste0(trait_name," (INT); 60-80% PGS groups; FOS"),x_name = 'PGS (scaled)',y_name = trait_name,width = 6, height = 5)

# plot scatter plot: 81-100% PGS; INT normalized phenotype; FOS_24h
df_to_plot <- subset(dt_s3_g5_INT, select=c("SCORE.norm", trait_field, "FOS_24h"))
colnames(df_to_plot) <- c("x", "y", "group")
df_to_plot$group[df_to_plot$group==1] <- "FOS"
df_to_plot$group[df_to_plot$group==0] <- "non-FOS"
plot_scatter1(png_name=paste0(outDir,"/scatter_strat_80To100PGS_INT_FOS_",trait_name,".png"),data = df_to_plot,title = paste0(trait_name," (INT); 80-100% PGS groups; FOS"),x_name = 'PGS (scaled)',y_name = trait_name,width = 6, height = 5)

####################
# plot boxplot: PGS percentiles; raw phenotype; Fish oil
df_to_plot <- rbind(subset(dt_s2_g1,select=c("SCORE.norm", trait_field, "FOS_24h", "pgs_perct")),
		    subset(dt_s2_g2,select=c("SCORE.norm", trait_field, "FOS_24h", "pgs_perct")),
		    subset(dt_s3_g1, select=c("SCORE.norm", trait_field, "FOS_24h", "pgs_perct")),
		    subset(dt_s3_g2, select=c("SCORE.norm", trait_field, "FOS_24h", "pgs_perct")),
		    subset(dt_s3_g3, select=c("SCORE.norm", trait_field, "FOS_24h", "pgs_perct")),
		    subset(dt_s3_g4, select=c("SCORE.norm", trait_field, "FOS_24h", "pgs_perct")),
		    subset(dt_s3_g5, select=c("SCORE.norm", trait_field, "FOS_24h", "pgs_perct")))
plot_dt_s1 <- copy(dt_s1); plot_dt_s1$pgs_perct <- "All"
df_to_plot <- rbind(df_to_plot, subset(plot_dt_s1, select=c("SCORE.norm", trait_field, "FOS_24h", "pgs_perct")))
df_to_plot$pgs_perct <- factor(df_to_plot$pgs_perct, levels=c("0-5%","0-20%","20-40%","40-60%","60-80%","80-100%","95-100%","All"))
colnames(df_to_plot) <- c("x", "y", "fo_group", "pgs_group")
df_to_plot$fo_group[df_to_plot$fo_group==1] <- "FOS"
df_to_plot$fo_group[df_to_plot$fo_group==0] <- "non-FOS"
plot_boxplot(png_name=paste0(outDir,"/boxplot_strat_PGSgroups_FOS_",trait_name,".png"),data = df_to_plot,title = paste0(trait_name," (raw); all PGS groups; FOS"), y_name = trait_name,width = 8, height = 6)

# plot boxplot: PGS percentiles; INT phenotype; Fish oil
df_to_plot <- rbind(subset(dt_s2_g1_INT,select=c("SCORE.norm", trait_field, "FOS_24h", "pgs_perct")),
                    subset(dt_s2_g2_INT,select=c("SCORE.norm", trait_field, "FOS_24h", "pgs_perct")),
                    subset(dt_s3_g1_INT, select=c("SCORE.norm", trait_field, "FOS_24h", "pgs_perct")),
                    subset(dt_s3_g2_INT, select=c("SCORE.norm", trait_field, "FOS_24h", "pgs_perct")),
                    subset(dt_s3_g3_INT, select=c("SCORE.norm", trait_field, "FOS_24h", "pgs_perct")),
                    subset(dt_s3_g4_INT, select=c("SCORE.norm", trait_field, "FOS_24h", "pgs_perct")),
                    subset(dt_s3_g5_INT, select=c("SCORE.norm", trait_field, "FOS_24h", "pgs_perct")))
plot_dt_s1_INT <- copy(dt_s1_INT); plot_dt_s1_INT$pgs_perct <- "All"
df_to_plot <- rbind(df_to_plot, subset(plot_dt_s1_INT, select=c("SCORE.norm", trait_field, "FOS_24h", "pgs_perct")))
df_to_plot$pgs_perct <- factor(df_to_plot$pgs_perct, levels=c("0-5%","0-20%","20-40%","40-60%","60-80%","80-100%","95-100%","All"))
colnames(df_to_plot) <- c("x", "y", "fo_group", "pgs_group")
df_to_plot$fo_group[df_to_plot$fo_group==1] <- "FOS"
df_to_plot$fo_group[df_to_plot$fo_group==0] <- "non-FOS"
plot_boxplot(png_name=paste0(outDir,"/boxplot_strat_PGSgroups_INT_FOS_",trait_name,".png"),data = df_to_plot,title = paste0(trait_name," (INT); all PGS groups; FOS"), y_name = trait_name,width = 8, height = 6)

}
