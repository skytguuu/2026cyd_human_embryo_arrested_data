# 6. M-decay基因分析

library(dplyr)
library(ggpubr)
library(tidyverse)

d <- read.table("GSE36552_nsb_hum_embry_genefpkm.txt", header=TRUE, sep="\t",
                quote="", comment.char="", fill=TRUE, check.names=FALSE)
gv = read.table("GSE107746_Folliculogenesis_FPKM.log2.txt", header=TRUE, stringsAsFactors=FALSE)

y1 = data.frame(geneid=d$Gene_ID,
                zyg=apply(d[, c(3:5)], 1, mean),
                cell8=apply(d[, c(6:dim(d)[2])], 1, mean))
y2 = data.frame(geneid=gv$gene,
                gv=apply(gv[, grep("Antral_follicle", x=colnames(gv))], 1, mean))
dat = merge(y1, y2, by="geneid")

df <- dat %>% filter(gv > 2) %>%
     mutate(across(c(gv, zyg, cell8), ~ log2(.x + 1)))

gene_clusters <- df %>% mutate(Cluster=case_when(
    gv > zyg + 1 & zyg <= cell8 + 1 ~ "Cluster I",
    gv <= zyg + 1 & gv > zyg - 1 & zyg > cell8 + 1 ~ "Cluster II",
    gv > zyg + 1 & zyg > cell8 + 1 ~ "Cluster III",
    gv <= zyg + 1 & gv > zyg - 1 & zyg <= cell8 + 1 & zyg > cell8 - 1 ~ "Cluster IV",
    TRUE ~ "Unclassified"))

gene_clusters$type = ifelse(gene_clusters$Cluster %in% c("Cluster II", "Cluster III"), "Z-decay",
                     ifelse(gene_clusters$Cluster == "Cluster I", "M-decay", gene_clusters$Cluster))
zdecay = subset(gene_clusters, gene_clusters$type=="Z-decay")
mdecay = subset(gene_clusters, gene_clusters$type=="M-decay")

zga = read.table("hum_8vs4_ZGA_fc2_p0.01.txt", header=TRUE)
dat = read.csv("hum_8cell_rnaseq_fpkm_symbol.csv", header=TRUE, stringsAsFactors=FALSE)

# 基因比较
plot_dat <- dat %>% filter(group %in% c("ZGA", "Z-decay", "M-decay")) %>%
  pivot_longer(cols=c("ab", "normal"), names_to="Status", values_to="Expression") %>%
  mutate(Expression=log2(Expression + 1)) %>%
  mutate(Status=factor(Status, levels=c("normal", "ab")),
         group=factor(group, levels=c("M-decay", "ZGA", "Z-decay")))

p <- ggplot(plot_dat, aes(x=Status, y=Expression, color=Status)) +
  geom_boxplot(outlier.shape=NA, width=0.6, lwd=0.8) +
  facet_wrap(~group, scales="free") +
  scale_color_manual(values=c("normal"="#74add1", "ab"="#f46d43")) +
  stat_compare_means(method="t.test", label="p.format", label.x=1.5) +
  theme_bw() +
  theme(
    strip.background=element_blank(),
    strip.text=element_text(size=12, face="bold"),
    axis.text.x=element_text(angle=45, hjust=1),
    legend.position="none",
    panel.grid=element_blank()
  ) +
  labs(x=NULL, y=expression(log[2](FPKM + 1)))
