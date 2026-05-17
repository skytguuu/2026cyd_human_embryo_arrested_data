# 1. 基本HiC数据分析 - R部分: 跨染色体散点图与比例柱状图

library(ggplot2)
library(reshape2)

ref = read.table("~/data1/Annotation/HiC_needfile/hicpro/chrom_hg19.sizes",
                 header=FALSE, sep="\t", stringsAsFactors=FALSE)[1:22,]
out = as.data.frame(t(combn(chrOrder, 2)))
tmp = cbind(chrOrder, chrOrder)
colnames(tmp) = colnames(out)
out = rbind(out, tmp)

for (i in 1:nrow(out)){
  out[i,3] = ref[ref$V1 %in% out[i,1], 2] / ref[ref$V1 %in% out[i,2], 2]
  out[i,4] = con_data[rownames(con_data)==out[i,1], colnames(con_data)==out[i,2]]
  out[i,5] = dif_data[rownames(dif_data)==out[i,1], colnames(dif_data)==out[i,2]]
}
colnames(out)[3:5] = c("ratio", "con_value", "dif_value")
out$ratio = log2(out$ratio)

# 散点拟合
p = ggplot(data=out[,-c(1,2)], aes(x=ratio, y=con_value), color="blue") +
    geom_point(size=3) + stat_smooth(method='lm')
label = cor(out$ratio, out$con_value)
label = paste0("Pearson's r=", label)
p = p + theme_bw() + geom_text(aes(x=1, y=2.5, label=label))
p

# 两组相比
out$dif2con = out$dif_value / out$con_value
p = ggplot(data=out[,-c(1,2)], aes(x=ratio, y=dif2con), color="blue") +
    geom_point(size=3) + stat_smooth(method='lm')
label = cor(out$ratio, out$dif2con)
label = paste0("Pearson's r=", label)
p = p + theme_bw() + geom_text(aes(x=1, y=2, label=label))
p
out$fc = out$dif_value / out$con_value

# 比例柱状图
con_intra = sum(diag(as.matrix(con_convert)))
dif_intra = sum(diag(as.matrix(dif_convert)))
res = data.frame(NC=c(con_intra, (sum(con_convert)-con_intra)/2),
                 SH=c(dif_intra, (sum(dif_convert)-dif_intra)/2))
res = apply(res, 2, function(x){x/sum(x)})
rownames(res) = c("intra", "inter")

mm = reshape2::melt(res)
ggplot(mm, aes(x=Var2, y=value, fill=factor(Var1, levels=c("inter", "intra")))) +
  geom_bar(position="fill", stat="identity") +
  geom_text(aes(label=scales::percent(value, accuracy=0.1)),
            position=position_fill(vjust=0.5)) +
  theme_bw() +
  scale_fill_manual(values=ggsci::pal_nejm()(2))
