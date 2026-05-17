# 5. Loop分析 - R部分

library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(stringr)

# 读取loop数据
nor = read.table("output1.txt", header=FALSE)
abnor = read.table("output2.txt", header=FALSE)

nor$symbol = stringr::str_split(nor$V4, pattern=":", simplify=TRUE)[, 1]
abnor$symbol = stringr::str_split(abnor$V4, pattern=":", simplify=TRUE)[, 1]

# GO富集
ID1 = bitr(nor$Symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")[, 2]
ego1 = enrichGO(gene=ID1, 'org.Hs.eg.db', ont="BP", pvalueCutoff=0.05,
                pAdjustMethod="BH", qvalueCutoff=0.05, maxGSSize=10000, readable=TRUE)
c1 = ego1@result[ego1@result$qvalue < 0.05,]
c1 = as.data.frame(c1)

# 差异染色质环交集差异基因
deg = read.table("../rnaseq/cyd_h8cell_nor2ab_limma_DEG.txt", header=TRUE)
nor.spe.deg = intersect(nor.spe, deg$Symbol)
abnor.spe.deg = intersect(abnor.spe, deg$Symbol)

ID1 = bitr(nor.spe.deg, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")[, 2]
ID2 = bitr(abnor.spe.deg, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")[, 2]
ego1 = enrichGO(gene=ID1, 'org.Hs.eg.db', ont="BP", pvalueCutoff=0.05,
                pAdjustMethod="BH", qvalueCutoff=0.05, maxGSSize=10000, readable=TRUE)
ego2 = enrichGO(gene=ID2, 'org.Hs.eg.db', ont="BP", pvalueCutoff=0.05,
                pAdjustMethod="BH", qvalueCutoff=0.05, maxGSSize=10000, readable=TRUE)

# 这些基因富集在哪个TAD中
nor.tad = read.table("../tad/normal_TAD_40kb_genes.txt", header=FALSE)
abnor.tad = read.table("../tad/abnormal_TAD_40kb_genes.txt", header=FALSE)
colnames(nor.tad) = c("seqnames", "start", "end", "genechr", "genestart", "geneend", "symbol")
colnames(abnor.tad) = c("seqnames", "start", "end", "genechr", "genestart", "geneend", "symbol")
nor.tad$id = paste(nor.tad$seqnames, nor.tad$start, nor.tad$end, sep="_")
abnor.tad$id = paste(abnor.tad$seqnames, abnor.tad$start, abnor.tad$end, sep="_")
nor.tad.ref = nor.tad[, c(8, 7)]
abnor.tad.ref = abnor.tad[, c(8, 7)]
out1 <- enricher(nor.spe.deg, TERM2GENE=nor.tad.ref)
out2 <- enricher(abnor.spe.deg, TERM2GENE=abnor.tad.ref)
g1 = unlist(stringr::str_split(out1@result$geneID, pattern="/"))
g2 = unlist(stringr::str_split(out2@result$geneID, pattern="/"))

# 画染色体loop环个数柱状图
nor.num <- nor %>%
    na.omit() %>%
    group_by(V1) %>%
    summarize(count=n())
nor.num$type = "normal"
inp = rbind(nor.num, abnor.num)
inp$V1 = factor(inp$V1, levels=c(paste0("chr", 1:22), "chrX"))

ggplot(inp, aes(x=V1, y=count, fill=type)) +
    geom_col(position="dodge") +
    labs(x="V1", y="Count") +
    scale_fill_manual(values=ggsci::pal_nejm()(2)) +
    theme_minimal() + theme_bw()
