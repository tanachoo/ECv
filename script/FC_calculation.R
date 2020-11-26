# Author: yoshi
# Project: GSE49644, PSB2020
# Write Start: 3/18/2019
# Updated: 7/16/2019

## script for calcucating log2FC in GSE49644 microarray data

# input data is log2-transformed
data <- read.table("GSE49644.txt", header=T, row.name=1, sep="\t", quote="", as.is=T, check.names=F)

colnames(data)

# change colnames
new_col <- c("HCC827_ctl_1", "HCC827_ctl_2", "HCC827_ctl_3", "A549_ctl_1", "A549_ctl_2", "A549_ctl_3", "NCIH358_ctl_1", "NCIH358_ctl_2", "NCIH358_ctl_3", "A549_EMT_1", "A549_EMT_2", "A549_EMT_3", "HCC827_EMT_1", "HCC827_EMT_2", "HCC827_EMT_3", "NCIH358_EMT_1", "NCIH358_EMT_2", "NCIH358_EMT_3")
colnames(data) <- new_col


# Add mean columns
#data_add <- within(data, HCC827_control_mean <- (data[, 1] + data[, 2] + data[, 3])/ 3)
#data_add <- within(data_add, A549_control_mean <- (data[, 4] + data[, 5] + data[, 6])/ 3)
#data_add <- within(data_add, NICH358_control_mean <- (data[, 7] + data[, 8] + data[, 9])/ 3)

# updated codes
data_add <- within(data, HCC827_control_mean <- (data[, "HCC827_ctl_1"] + data[, "HCC827_ctl_2"] + data[, "HCC827_ctl_3"])/ 3)
data_add <- within(data_add, A549_control_mean <- (data[, "A549_ctl_1"] + data[, "A549_ctl_2"] + data[, "A549_ctl_3"])/ 3)
data_add <- within(data_add, NICH358_control_mean <- (data[, "NCIH358_ctl_1"] + data[, "NCIH358_ctl_2"] + data[, "NCIH358_ctl_3"])/ 3)


#data_add <- within(data_add, HCC827_EMT_mean <- (data[, 13] + data[, 14] + data[, 15])/ 3)
#data_add <- within(data_add, A549_EMT_mean <- (data[, 10] + data[, 11] + data[, 12])/ 3)
#data_add <- within(data_add, NCIH358_EMT_mean <- (data[, 16] + data[, 17] + data[, 18])/ 3)

# updated codes
data_add <- within(data_add, HCC827_EMT_mean <- (data[, "HCC827_EMT_1"] + data[, "HCC827_EMT_2"] + data[, "HCC827_EMT_3"])/ 3)
data_add <- within(data_add, A549_EMT_mean <- (data[, "A549_EMT_1"] + data[, "A549_EMT_2"] + data[, "A549_EMT_3"])/ 3)
data_add <- within(data_add, NCIH358_EMT_mean <- (data[, "NCIH358_EMT_1"] + data[, "NCIH358_EMT_2"] + data[, "NCIH358_EMT_3"])/ 3)

write.table(data_add, "GSE49644_add_mean.txt", sep="\t", append=F, quote=F, row.names=T, col.names=NA)

# Add log2FC columns
#data_add <- within(data_add, HCC827_FC <- (data_add[, 22] - data_add[, 19]))
#data_add <- within(data_add, A549_FC <- (data_add[, 23] - data_add[, 20]))
#data_add <- within(data_add, NCIH358_FC <- (data_add[, 24] - data_add[, 21]))

# updated codes
data_add <- within(data_add, HCC827_FC <- (data_add[, "HCC827_EMT_mean"] - data_add[, "HCC827_control_mean"]))
data_add <- within(data_add, A549_FC <- (data_add[, "A549_EMT_mean"] - data_add[, "A549_control_mean"]))
data_add <- within(data_add, NCIH358_FC <- (data_add[, "NCIH358_EMT_mean"] - data_add[, "NICH358_control_mean"]))

write.table(data_add, "GSE49644_add.txt", sep="\t", append=F, quote=F, row.names=T, col.names=NA)


