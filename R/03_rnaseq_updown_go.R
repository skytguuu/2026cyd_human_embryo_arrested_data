# 3. RNAseq分析 - 上下调基因GO与GSEA

library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# 上下调基因
up = subset(res_1, res_1$logFC > 0)
down = subset(res_1, res_1$logFC < 0)
ID1 = bitr(up$Symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")[, 2]
ego1 = enrichGO(gene=ID1, 'org.Hs.eg.db', ont="BP", pvalueCutoff=0.05,
                pAdjustMethod="BH", qvalueCutoff=0.05, maxGSSize=10000, readable=TRUE)
ID2 = bitr(down$Symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")[, 2]
ego2 = enrichGO(gene=ID2, 'org.Hs.eg.db', ont="BP", pvalueCutoff=0.05,
                pAdjustMethod="BH", qvalueCutoff=0.05, maxGSSize=10000, readable=TRUE)

c1 = ego1@result[ego1@result$qvalue < 0.05,]
c1 = as.data.frame(c1)
path = rbind(as.data.frame(ego1@result), as.data.frame(ego2@result)) %>%
       select(c("ID", "p.adjust", "qvalue", "Count"))
path$source = c(rep("up", dim(ego1@result)[1]), rep("down", dim(ego2@result)[1]))
path = subset(path, path$p.adjust < 0.05)
path$fdr = log10(path$p.adjust)
path$fdr = ifelse(path$source=="down", path$fdr, -path$fdr)

ggplot(path, aes(x=reorder(ID, fdr, decreasing=FALSE), y=fdr, fill=source)) +
    geom_col(width=0.7) + coord_flip() +
    scale_fill_manual(values=c("up"="#FFCCCC", "down"="#CCCCFF")) +
    theme_minimal() +
    theme(
        axis.title.y=element_blank(),
        axis.text=element_text(color="black", face="bold"),
        panel.grid=element_blank(),
        axis.line.x=element_line(color="black", size=0.5),
        axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position="right",
        plot.title=element_text(hjust=0.5, face="bold")
    ) + labs(x=NULL, y="-log10 (FDR)", title="C7") +
      scale_y_continuous(limits=c(-8, 16), breaks=seq(-8, 16, 4))

# GSEA
g1 = res_1[order(res_1$logFC, decreasing=TRUE), ]
geneList = g1$logFC
names(geneList) = g1$Symbol
names(geneList) = bitr(names(geneList), fromType="SYMBOL", toType="ENTREZID",
                       OrgDb="org.Hs.eg.db")[, 2]

ego3 <- gseGO(geneList=geneList,
              OrgDb=org.Hs.eg.db,
              ont="BP",
              minGSSize=100,
              maxGSSize=20000,
              pvalueCutoff=0.05,
              verbose=FALSE)

library(msigdbr)
m_t2g <- msigdbr(species="Homo sapiens", category="C2", subcategory="CGP") %>%
         dplyr::select(gs_name, gene_symbol)
em2 <- GSEA(geneList, TERM2GENE=m_t2g)

library(enrichplot)
library(RColorBrewer)
source("/home/dell/data1/softwares/Rscript/gseaplot4.R")
gseaplot4(em2, geneSetID=1:3, pvalue_table=TRUE,
          color=c("#E495A5", "#86B875", "#7DB0DD"), ES_geom="dot")
