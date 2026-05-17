# 3. RNAseq分析 - 差异基因与通路富集

library(limma)
library(edgeR)
library(dplyr)
library(data.table)
library(org.Hs.eg.db)
options(stringsAsFactors=FALSE)

d = read.table("./hum_8cell_merge_RNAseq1.txt", header=TRUE, sep="\t")
input = d[, c(1, 7:8, 5:6)]
input$id <- gsub("\\.\\d+$", "", input$tracking_id)

# ID转化
input$symbol <- mapIds(org.Hs.eg.db, keys=as.character(input$id),
                        column="SYMBOL", keytype="ENSEMBL", multiVals="first")
r2 <- input[, c(2:5, 7)] %>%
  na.omit() %>%
  group_by(symbol) %>%
  summarize_all(mean)
r2 = as.data.frame(r2)
expr = r2[, -1]
rownames(expr) = r2[, 1]

# 差异基因
group <- factor(rep(c('control', 'treat'), each=2), levels=c('control', 'treat'))
design <- model.matrix(~0+group)
rownames(design) = colnames(expr)
colnames(design) <- levels(group)
fit <- lmFit(expr, design)
cont.matrix <- makeContrasts(contrasts=c('treat-control'), levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
nrDEG_limma_voom = topTable(fit2, coef='treat-control', n=Inf)
nrDEG_limma_voom = na.omit(nrDEG_limma_voom)
head(nrDEG_limma_voom)

library(dplyr)
res <- cbind(rownames(nrDEG_limma_voom), nrDEG_limma_voom)
res_1 <- res %>% dplyr::filter((logFC > log2(4) | logFC < (-log2(4))) & adj.P.Val < 0.05)
colnames(res_1)[1] <- "Symbol"
write.table(res_1, "cyd_h8cell_nor2ab_limma_DEG.txt", quote=FALSE, sep="\t", row.names=FALSE)

# 差异通路
library(clusterProfiler)
ID1 = bitr(res_1$Symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")[, 2]
ego1 = enrichGO(gene=ID1, 'org.Hs.eg.db', ont="BP", pvalueCutoff=0.05,
                pAdjustMethod="BH", qvalueCutoff=0.05, maxGSSize=10000, readable=TRUE)
c1 = ego1@result[ego1@result$qvalue < 0.05,]
c1 = as.data.frame(c1)
write.csv(c1, "cyd_deg_go.csv", quote=FALSE, row.names=FALSE)

d = read.csv("cyd_deg_go.csv", header=TRUE)
library(ggplot2)
d$Description = factor(d$Description, ordered=FALSE)
p1 = ggplot(d, aes(richfactor, Description))
p2 = p1 + geom_point()
pbubble = p2 + geom_point(aes(size=Count, color=-1*log10(qvalue))) + scale_size(range=c(4, 10))
pr = pbubble + scale_color_gradient2(low='#0100FF', mid='#6C02E2', high='#FF0101', midpoint=10)
pr = pr + labs(color=expression(-log[10](FDR)), size="Gene number",
               x="Rich factor", y="Top of 10 GO enrichment") +
          theme(plot.title=element_text(hjust=0.5))
pr = pr + theme_bw()
pr
