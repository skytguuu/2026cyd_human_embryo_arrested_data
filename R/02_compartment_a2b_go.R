# 2. AB区室 - R部分: 方法2 - AB转化区域GO富集分析

library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

color <- brewer.pal(3, "Dark2")
colorl <- rep(color, each=10)

a2b = read.table("h8cell_a2b.genes")
a2b$gene = sub("\\.\\d+$", "", a2b$V4)
ID1 = bitr(a2b$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")[, 2]
up1 <- enrichGO(gene=ID1, OrgDb="org.Hs.eg.db", ont="BP",
                pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.05,
                readable=TRUE, maxGSSize=10000)
dotplot(up1, showCategory=10)

c1 = up1@result[up1@result$qvalue < 0.05,]
c1$type = "A2B"
c1$source = "FGSC2GV"
path = rbind(c1, c2, m1, m2)
pdata = path %>% group_by(source, type) %>% top_n(n=-10, wt=p.adjust)
pdata1 = subset(pdata, pdata$source=='FGSC2GV')
pdata1$Description = factor(pdata1$Description, levels=pdata1$Description)

ggplot(pdata1, aes(x=Description, y=-log10(p.adjust), fill=type)) +
  geom_bar(stat="identity") + coord_flip() +
  scale_fill_manual(values=c("A2B"="#2E8B57", "B2A"="#D2691E")) +
  scale_color_manual(values=c("A2B"="#2E8B57", "B2A"="#D2691E")) +
  theme(
    axis.title=element_text(size=15, face="plain", color="black"),
    axis.text=element_text(size=12, face="plain", color="black"),
    axis.text.y=element_text(color=ifelse(pdata1$type=="A2B", "#2E8B57", "#D2691E")),
    axis.title.y=element_blank(),
    legend.title=element_blank(),
    legend.text=element_text(size=8, face="bold"),
    legend.margin=margin(t=0, r=0, b=0, l=0, unit="pt"),
    legend.direction="vertical",
    legend.position=c(0.9, 0.92),
    legend.background=element_blank(),
    panel.background=element_rect(fill="transparent", colour="black"),
    plot.background=element_blank()
  )
