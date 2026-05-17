# 2. AB区室 - R部分: 使用HiTC计算AB区室

library("HiTC")

name = list.files(pattern="_500000_iced.matrix")
n = length(name)

pas = function(x){
  a = paste(x[,1], x[,2], sep=":")
  b = paste(a, x[,3], sep="-")
  return(b)
}

M = c()
for (j in 1:22){
    chr.name = paste("chr", j, sep="")
    chr_matrix = c()
    for (i in 1:n){
        Exp = importC(name[i],
                      xgi.bed="~/data1/Annotation/HiC_needfile/human_genome_500kb.bed",
                      allPairwise=FALSE, rm.trans=FALSE, lazyload=FALSE)
        if (i == 1){
            chr.name <- pca.hic(Exp[isIntraChrom(Exp)][[j]],
                                normPerExpected=TRUE, method="mean", npc=1)
            chr.name = as.matrix(as.data.frame(chr.name))
            chr.name = chr.name[, c(3:5, 8)]
            colnames(chr.name)[4] = name[i]
            chr_matrix = chr.name
        } else {
            chr.name <- pca.hic(Exp[isIntraChrom(Exp)][[j]],
                                normPerExpected=TRUE, method="mean", npc=1)
            chr.name = as.matrix(as.data.frame(chr.name))
            chr.name = chr.name[, c(3:5, 8)]
            colnames(chr.name)[4] = name[i]
            chr_matrix = merge(chr_matrix, chr.name, by=c("seqnames", "start", "end"))
        }
    }
    M = rbind(M, chr_matrix)
}

write.table(M, "cyd_chr1-22_mean_500kb.txt", row.names=FALSE, quote=FALSE,
            col.names=TRUE, sep="\t")

library(stringr)
id = str_split(name, pattern="_merge", simplify=TRUE)[,1]
for (i in 4:ncol(M)){
    label = paste(id[i-3], "_chr1-22_mean_500kb.txt", sep="")
    write.table(M[, c(1:3, i)], label, row.names=FALSE, quote=FALSE,
                col.names=FALSE, sep="\t")
}
