# 6. Minor ZGA分析

library(limma)
library(edgeR)
library(dplyr)
library(org.Hs.eg.db)

d1 = read.table("GSE36552_nsb_hum_embry_genefpkm.txt", header=TRUE, sep="\t",
                quote="", comment.char="", fill=TRUE, check.names=FALSE)
dat = read.csv("hum_8cell_rnaseq_fpkm_symbol.csv", header=TRUE, stringsAsFactors=FALSE)

# 匹配包含 "2.cell" 或 "4.cell" 的列名
target_columns <- grep("Gene_ID|2.cell|4.cell", colnames(d1), value=TRUE)
d_subset <- d1[, target_columns]

input = d_subset[, c(1:19)]
expr = as.data.frame(input[, -1])
rownames(expr) = input$Gene_ID

# limma找差异
group <- factor(rep(c('control', 'treat'), times=c(6, 12)), levels=c('control', 'treat'))
design <- model.matrix(~0+group)
rownames(design) = colnames(expr)
colnames(design) <- levels(group)
fit <- lmFit(expr, design)
cont.matrix <- makeContrasts(contrasts=c('treat-control'), levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
nrDEG_limma_voom = topTable(fit2, coef='treat-control', n=Inf)
nrDEG_limma_voom = na.omit(nrDEG_limma_voom)

library(dplyr)
res <- cbind(rownames(nrDEG_limma_voom), nrDEG_limma_voom)
res_1 <- res %>% dplyr::filter((logFC > 2 | logFC < (-2)) & adj.P.Val < 0.01)
colnames(res_1)[1] <- "Symbol"
ZGA = subset(res_1, res_1$logFC > 0)
write.table(ZGA, "hum_4vs2_minorZGA_fc2_p0.0", quote=FALSE, sep="\t", row.names=FALSE)

# 画图
d3 = dat[dat$symbol %in% rownames(ZGA), ]
