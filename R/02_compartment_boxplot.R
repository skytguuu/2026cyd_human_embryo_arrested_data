# 2. AB区室 - R部分: Compartment strength箱线图

library("ggpubr")
library(reshape2)

d = read.table("compartment_saddle_strength.txt", header=TRUE, sep="\t")
d1 = as.data.frame(t(d[,-1]))
dd = melt(d1)

p = ggplot(dd, aes(x=dd$variable, y=dd$value, fill=dd$variable)) +
    geom_boxplot(position=position_dodge(0.8))
p = p + stat_compare_means(comparisons=combn(levels(dd$variable), 2, simplify=FALSE)) + theme_bw()
p = p + stat_compare_means(comparisons=c("Abnormal", "Normal")) + theme_bw()
vioplot::vioplot(d1, col=c("#E69F00", "#56B4E9", "#CC79A7", "#009E73"))
