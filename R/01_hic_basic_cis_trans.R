# 1. 基本HiC数据分析 - R部分: cis/trans比例箱线图

library(ggplot2)
library(ggpubr)

nc = read.table("hum_8cell_2PN_merge.allValidPairs.txt", header=FALSE, stringsAsFactors=FALSE)
sh = read.table("hum_8cell_abnorm_2PN_merge.allValidPairs.txt", header=FALSE, stringsAsFactors=FALSE)
boxplot(nc$V4, sh$V4, outline=FALSE, col="orange")
