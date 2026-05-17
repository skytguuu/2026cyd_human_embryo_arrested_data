# 4. TAD分析 - R部分: TAD RTI与边界分析

# TAD RTI
normal = read.csv("./normal_RTI_40kb_TAD.csv", header=TRUE, stringsAsFactors=FALSE)
abnormal = read.csv("./abnormal_RTI_40kb_TAD.csv", header=TRUE, stringsAsFactors=FALSE)
boxplot(normal$RTI_score, abnormal$RTI_score, outline=FALSE, col=rainbow(2))
vioplot::vioplot(normal$RTI_score, abnormal$RTI_score,
                 col=c("#E69F00", "#56B4E9"), ylab="Relative TAD score",
                 names=c("nor", "abn"))

# TAD boundary - 边界强度
normal = read.table("hum_8cell_2PN_40kb_input.is520001.ids400001.insulation.boundaries.bed",
                    header=FALSE, skip=1)
abnormal = read.table("hum_8cell_abnorm_2PN_40kb_input.is520001.ids400001.insulation.boundaries.bed",
                      header=FALSE, skip=1)
boxplot(normal$V5, abnormal$V5, outline=FALSE, col=rainbow(2),
        ylab="Boundary Score", xlab=c("Normal", "Abnormal"))

# TAD内部基因与表达
d = read.table("../compartment/500kb/hum_8cell_merge_RNAseq2.txt", header=TRUE, sep="\t")
nor.tad = read.table("./normal_TAD_40kb_genes.txt", header=FALSE)
abnor.tad = read.table("./abnormal_TAD_40kb_genes.txt", header=FALSE)
colnames(nor.tad) = c("seqnames", "start", "end", "genechr", "genestart", "geneend", "symbol")
colnames(abnor.tad) = c("seqnames", "start", "end", "genechr", "genestart", "geneend", "symbol")
com1 = merge(nor.tad, d, by="symbol")
com2 = merge(abnor.tad, d, by="symbol")

# 计算每个TAD的基因表达
library(dplyr)
nor.r2 <- com1[, c(2:4, 12:13)] %>%
    na.omit() %>%
    group_by(seqnames, start, end) %>%
    summarize_all(mean, na.rm=TRUE)
nor.r2 = as.data.frame(nor.r2)
abnor.r2 <- com2[, c(2:4, 12:13)] %>%
    na.omit() %>%
    group_by(seqnames, start, end) %>%
    summarize_all(mean, na.rm=TRUE)
abnor.r2 = as.data.frame(abnor.r2)
boxplot(nor.r2$normal, abnor.r2$abnor, outline=FALSE,
        ylab="TAD FPKM", names=c("normal", "abnormal"))
