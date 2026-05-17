# 2. AB区室 - R部分: AB区室基因与表达分析

library(org.Hs.eg.db)
library(dplyr)

d = read.table("../../rnaseq/hum_8cell_merge_RNAseq.txt", header=TRUE, sep="\t")
d$normal = apply(d[, 7:8], 1, mean)
d$abnor = apply(d[, 5:6], 1, mean)
d$id <- gsub("\\.\\d+$", "", d$tracking_id)

# ID转化
d$symbol <- mapIds(org.Hs.eg.db, keys=as.character(d$id),
                    column="SYMBOL", keytype="ENSEMBL", multiVals="first")
r2 <- d[, c(5:10, 12)] %>%
   na.omit() %>%
   group_by(symbol) %>%
   summarize_all(mean)
r2 = as.data.frame(r2)
write.table(r2, "hum_8cell_merge_RNAseq2.txt", row.names=FALSE, sep="\t")

# AB区室与基因表达合并
ref = read.table("gene_loated_500kb_compartment.txt", header=FALSE)
colnames(ref) = c("seqnames", "start", "end", "val", "genechr", "genestart", "geneend", "symbol")
compart = read.table("cyd_chr1-22_mean_500kb.txt", header=TRUE)
com1 = merge(compart, ref, by=colnames(ref)[1:3], all=TRUE)
com2 = merge(com1, r2, by="symbol")
com2$ratio = (com2$abnor + 0.01) / (com2$normal + 0.01)

normal.a = subset(com2, com2$hum_8cell_2PN_merge_500000_iced.matrix > 0)
normal.b = subset(com2, com2$hum_8cell_2PN_merge_500000_iced.matrix < 0)
boxplot(abnor.a$normal, abnor.b$normal, outline=FALSE, col="darkblue")
vioplot::vioplot(log2(normal.a$normal+0.1), log2(normal.b$normal+0.1),
                 log2(abnor.a$abnor+0.1), log2(abnor.b$abnor+0.1),
                 col=rep(c("darkblue", "#E69F00"), 2))

# A2B基因通路
com3 = com2
com3 = com2[com2$normal != 0,]
com3$ratio = com3$abnor / com3$normal
a2b = subset(com3, com3$hum_8cell_2PN_merge_500000_iced.matrix > 0 &
                 com3$hum_8cell_abnorm_2PN_merge_500000_iced.matrix < 0)
b2a = subset(com3, com3$hum_8cell_2PN_merge_500000_iced.matrix < 0 &
                 com3$hum_8cell_abnorm_2PN_merge_500000_iced.matrix > 0)
stable = subset(com3, com3$hum_8cell_2PN_merge_500000_iced.matrix *
                    com3$hum_8cell_abnorm_2PN_merge_500000_iced.matrix > 0)
boxplot(a2b$ratio, b2a$ratio, stable$ratio, outline=FALSE)

# A2B GO富集
ID1 = bitr(a2b$symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")[, 2]
a2b.go <- enrichGO(gene=ID1, OrgDb="org.Hs.eg.db", ont="BP",
                    pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.05,
                    readable=TRUE, maxGSSize=10000)
dotplot(a2b.go, showCategory=20)
c1 = a2b.go@result[a2b.go@result$qvalue < 0.05,]
c1[grep("chro", c1$Description),]
