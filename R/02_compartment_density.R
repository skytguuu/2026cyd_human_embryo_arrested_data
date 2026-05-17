# 2. AB区室 - R部分: PC密度图与散点密度图

# 画scatter密度图
library(MASS)
library(ggplot2)

d = read.table("cyd_chr1-22_mean_500kb.txt", header=TRUE, stringsAsFactors=FALSE)
d[is.na(d)] = 0
data = na.omit(d)

ggplot(data, aes(x=hum_8cell_2PN_merge_500000_iced.matrix,
                 y=hum_8cell_abnorm_2PN_merge_500000_iced.matrix)) +
  stat_density_2d(aes(fill=..level.., alpha=..level..), geom="polygon") +
  geom_point(size=0.6, alpha=0.3, color="lightgrey") +
  scale_fill_distiller(palette="Spectral") + theme_bw()

# density plot
library("ggplot2")
library("plyr")
library("reshape2")

d = read.table("cyd_chr1-22_mean_500kb.txt", header=TRUE, stringsAsFactors=FALSE)
d[is.na(d)] = 0
colnames(d)[4:5] = c("Normal", "Abnormal")
m = melt(d[, 4:5])

# statistic
library(dplyr)
stat = group_by(m, variable) %>% summarize(mean=median(value))

res = ggplot(m, aes(x=value, colour=variable)) +
    geom_density(aes(fill=variable), alpha=0.4, size=1.1) +
    geom_hline(yintercept=0, color="white", size=1.1)
res = res + theme_bw() + theme(panel.grid=element_blank(),
                                axis.line=element_line(colour="black"))
res1 = res + scale_fill_manual(values=c("#5F9BC6", "#FE9B56")) +
             scale_color_manual(values=c("#5F9BC6", "#FE9B56"))
res2 = res1 + geom_vline(data=stat, aes(xintercept=mean, color=variable),
                          linetype="dashed", size=1)
res2 = res2 + annotate("text", x=0.8, y=3, colour="#5F9BC6",
             label=as.character(paste("median", round(stat$mean[1], 2), sep="="))) +
             annotate("text", x=0.75, y=3.5, colour="#FE9B56",
             label=as.character(paste("median", round(stat$mean[2], 2), sep="=")))
res2
