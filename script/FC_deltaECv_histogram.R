
# Author: yoshi
# Updated: 7/30/2019
# Project: GSE49644, PSB2020

## Script for making histogram for ΔECv and FC
## NOTE: You need to change cell line names for A549, NCIH358, HCC827

library(MASS)
setwd("/Users/yoshi/Desktop/microarray_data/GSE49644/data")

ECv_data <- read.table("../BayesianNetwork_analysis/ECv/ECv_GSE49644_th010_list_merge.txt", header=T, row.name=NULL, sep="\t", quote="", as.is=T, check.names=F)

FC_data <- read.table("GSE49644_log2FC.txt", header=T, row.name=1, sep="\t", quote="", as.is=T, check.names=F)

'''
## Pattern 1
par(mar = c(4.5, 4.5, 4.5, 5.5))

truehist(abs(FC_data[, "HCC827_FC"]), xlim=c(0,10), prob=T,
         col="#ff990080", border="#ff990080", axes=FALSE, xlab="", ylab="")
axis(side=1)
axis(side = 2, col.axis = "#ff9900", col = "#ff9900")
mtext("log2FC", side = 2, line = 3)

par(new = TRUE)

truehist(abs(ECv_data[, "ECv_delta_abs_HCC827_EMT_ctl"]), xlim=c(0,10), prob=T,
         col="#66990080", border="#66990080", axes=FALSE, xlab="", ylab="")
axis(side = 4, col.axis = "#669900", col = "#669900")
mtext("ΔECv", side = 4, line = 3)
'''

## Pattern 2
breaks_ECv = seq(0, 5.0, 0.1) #HCC827(0, 5.5), A549(0, 6.5), NCIH358(0, 5.0)
breaks_FC = seq(0, 7.5, 0.1) #HCC827(0, 9), A549(0, 10.5), NCIH358(0, 7.5)

png("histogram_NCIH358.png", height = 1500, width = 1500, res = 200)
# green # need to change->xlim, xaxp
hist(abs(ECv_data[, "ECv_delta_abs_NCIH358_EMT_ctl"]), xlim=c(0, 8),　ylim=c(0, 8), xaxp=c(0, 8, 8),
     col="#00cc3350", border="#00cc33", xaxs="i", yaxs="i", las=1,
     xlab="threshold", ylab="density", main="NCI-H358", cex.main=1.5, cex.lab=1.4, cex.axis=1.3,
     breaks=breaks_ECv, freq=F)
# magenta
hist(abs(FC_data[, "NCIH358_FC"]),
     col="#ff00ff50", border="#ff00ff",
     breaks=breaks_FC, freq=F, add=T)

legend("topright", legend=c("ΔECv", "log2FC"), col=c("#00cc33", "#ff00ff"), pch=15)

dev.off()
