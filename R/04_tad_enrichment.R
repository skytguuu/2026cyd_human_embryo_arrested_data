# 4. TAD分析 - R部分: TAD富集分析

library(clusterProfiler)
library(dplyr)

deg = read.table("../rnaseq/cyd_h8cell_nor2ab_limma_DEG.txt", header=TRUE)
nor.tad = read.table("./normal_TAD_40kb_genes.txt", header=FALSE)
abnor.tad = read.table("./abnormal_TAD_40kb_genes.txt", header=FALSE)
colnames(nor.tad) = c("seqnames", "start", "end", "genechr", "genestart", "geneend", "symbol")
colnames(abnor.tad) = c("seqnames", "start", "end", "genechr", "genestart", "geneend", "symbol")
nor.tad$id = paste(nor.tad$seqnames, nor.tad$start, nor.tad$end, sep="_")
abnor.tad$id = paste(abnor.tad$seqnames, abnor.tad$start, abnor.tad$end, sep="_")
nor.tad.ref = nor.tad[, c(8, 7)]
abnor.tad.ref = abnor.tad[, c(8, 7)]

# 超几何检验
out1 <- enricher(deg$Symbol, TERM2GENE=abnor.tad.ref)

# 有基因表达差异的TAD区域的基因GO富集
gg = nor.tad.ref[grep("chr2_196960000_198640000", nor.tad.ref$id), ]$symbol
ID1 = bitr(gg, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")[, 2]
ego1 = enrichGO(gene=ID1, 'org.Hs.eg.db', ont="BP", pvalueCutoff=0.05,
                pAdjustMethod="BH", qvalueCutoff=0.05, maxGSSize=10000, readable=TRUE)

# GSEA方法
res = read.table("../rnaseq/cyd_h8cell_nor2ab_limma_all.txt", header=TRUE)
colnames(res)[1] = "Symbol"
g1 = res[order(res$logFC, decreasing=TRUE), ]
geneList = g1$logFC
names(geneList) = g1$Symbol
em2 <- GSEA(geneList, TERM2GENE=nor.tad.ref)
