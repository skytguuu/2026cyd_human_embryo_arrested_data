# 5. Loop分析 - R部分: AB区室转化的loop基因

ref = read.table("../compartment/500kb/gene_loated_500kb_compartment.txt", header=FALSE)
colnames(ref) = c("seqnames", "start", "end", "val", "genechr", "genestart", "geneend", "symbol")
compart = read.table("../compartment/500kb/cyd_chr1-22_mean_500kb.txt", header=TRUE)
com1 = merge(compart, ref, by=colnames(ref)[1:3])

a2b = subset(com1, com1$hum_8cell_2PN_merge_500000_iced.matrix > 0 &
                 com1$hum_8cell_abnorm_2PN_merge_500000_iced.matrix < 0)
b2a = subset(com1, com1$hum_8cell_2PN_merge_500000_iced.matrix < 0 &
                 com1$hum_8cell_abnorm_2PN_merge_500000_iced.matrix > 0)
stable = subset(com1, com1$hum_8cell_2PN_merge_500000_iced.matrix *
                    com1$hum_8cell_abnorm_2PN_merge_500000_iced.matrix > 0)

# loop基因
nor.spe.a2b = intersect(a2b$symbol, nor.spe)
abnor.spe.b2a = intersect(b2a$symbol, abnor.spe)
