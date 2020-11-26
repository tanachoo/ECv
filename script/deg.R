# Author: yoshi
# Date: 12/13/2019
# Updated: 
# Project: ECv
# Script: Extract Differential Express Genes for GSE49644 dataset
# Extract DEG for each cell lines (A549, HCC827, NCI-H358) under condition between ctl VS TGFβ


## Load data
data <- read.table("./GSE49644.txt", header=T, 
                    row.name=1, sep="\t", quote="", as.is=T, check.names=F)

## Prep column for each cell line
col_HCC827 <- c('638297.HCC827.Parental', '638306.HCC827.Parental', '638307.HCC827.Parental',
                '638308.HCC827.EMT', '638309.HCC827.EMT', '638310.HCC827.EMT')

col_A549 <- c('638298.A549.Parental', '638301.A549.Parental', '638302.A549.Parental',
              '638303.A549.EMT', '638304.A549.EMT', '638305.A549.EMT')

col_NCIH358 <- c('638299.NCI-H358.Parental', '638311.NCI-H358.Parental', '638312.NCI-H358.Parental',
                 '638313.NCI-H358.EMT', '638314.NCI-H358.EMT', '638315.NCI-H358.EMT')

## Cut data for each cell line
data_HCC827 <- data[, col_HCC827] # genes(19849) × samples(6)
data_A549 <- data[, col_A549]
data_NCIH358 <- data[, col_NCIH358]


## DEG for HCC827
library('limma')
design <- cbind(CTL=1, CTLvsTreat=c(0,0,0,1,1,1)) # design matrix
fit <- lmFit(data_HCC827, design)
fit <- eBayes(fit)
#topTable(fit, coef='CTLvsTreat', adjust='BH')
p_value <- fit$p.value[, 2] #fitのp.valueの2列目にp値
q_value <- p.adjust(p_value, method='BH')
ranking <- rank(p_value)

CTLmean_HCC827 <- (data_HCC827[, '638297.HCC827.Parental']
                   + data_HCC827[, '638306.HCC827.Parental']
                   + data_HCC827[, '638307.HCC827.Parental'])/3

EMTmean_HCC827 <- (data_HCC827[, '638308.HCC827.EMT']
                   + data_HCC827[, '638309.HCC827.EMT']
                   + data_HCC827[, '638310.HCC827.EMT'])/3

log2FC <- EMTmean_HCC827 - CTLmean_HCC827

output <- cbind(rownames(data_HCC827), data_HCC827, log2FC, p_value, q_value, ranking)
write.table(output, './GSE49644_DEG_HCC827.txt', sep='\t', append=F, quote=F, row.names=F)



## DEG for A549
design <- cbind(CTL=1, CTLvsTreat=c(0,0,0,1,1,1)) # design matrix
fit <- lmFit(data_A549, design)
fit <- eBayes(fit)
#topTable(fit, coef='CTLvsTreat', adjust='BH')
p_value <- fit$p.value[, 2] #fitのp.valueの2列目にp値
q_value <- p.adjust(p_value, method='BH')
ranking <- rank(p_value)

CTLmean_A549 <- (data_A549[, '638298.A549.Parental']
                 + data_A549[, '638301.A549.Parental']
                 + data_A549[, '638302.A549.Parental'])/3

EMTmean_A549 <- (data_A549[, '638303.A549.EMT']
                 + data_A549[, '638304.A549.EMT']
                 + data_A549[, '638305.A549.EMT'])/3

log2FC <- EMTmean_A549 - CTLmean_A549

output <- cbind(rownames(data_A549), data_A549, log2FC, p_value, q_value, ranking)
write.table(output, './GSE49644_DEG_A549.txt', sep='\t', append=F, quote=F, row.names=F)


## DEG for NCI-H358
design <- cbind(CTL=1, CTLvsTreat=c(0,0,0,1,1,1)) # design matrix
fit <- lmFit(data_NCIH358, design)
fit <- eBayes(fit)
#topTable(fit, coef='CTLvsTreat', adjust='BH')
p_value <- fit$p.value[, 2] #fitのp.valueの2列目にp値
q_value <- p.adjust(p_value, method='BH')
ranking <- rank(p_value)

CTLmean_NCIH358 <- (data_NCIH358[, '638299.NCI-H358.Parental']
                    + data_NCIH358[, '638311.NCI-H358.Parental']
                    + data_NCIH358[, '638312.NCI-H358.Parental'])/3

EMTmean_NCIH358 <- (data_NCIH358[, '638313.NCI-H358.EMT']
                   + data_NCIH358[, '638314.NCI-H358.EMT']
                   + data_NCIH358[, '638315.NCI-H358.EMT'])/3

log2FC <- EMTmean_NCIH358 - CTLmean_NCIH358

output <- cbind(rownames(data_NCIH358), data_NCIH358, log2FC, p_value, q_value, ranking)
write.table(output, './GSE49644_DEG_NCIH358.txt', sep='\t', append=F, quote=F, row.names=F)




