# 1. 基本HiC数据分析 - R部分: 跨染色体互作热图

library(pheatmap)
library(reshape2)
library(data.table)
library(ggplot2)

con = read.table("hum_8cell_2PN_merge_inter_chr.txt", header=FALSE, stringsAsFactors=FALSE)
dif = read.table("hum_8cell_abnorm_2PN_merge_inter_chr.txt", header=FALSE, stringsAsFactors=FALSE)

# 长数据转化成宽数据
chrOrder = paste("chr", 1:22, sep="")
con_convert = dcast(con, factor(V1, levels=chrOrder)~factor(V2, levels=chrOrder))[,-1]
dif_convert = dcast(dif, factor(V1, levels=chrOrder)~factor(V2, levels=chrOrder))[,-1]
rownames(con_convert) = rownames(dif_convert) = chrOrder

# 补齐对称矩阵
con_convert[lower.tri(con_convert)] = t(con_convert)[lower.tri(con_convert)]
dif_convert[lower.tri(dif_convert)] = t(dif_convert)[lower.tri(dif_convert)]

# 计算Exp矩阵
# 根据公式sum(diag(j))/number of diag(j)
obs_exp_cal = function(file){
  n = nrow(file)
  for (i in 0:(n-1)){
    label = (row(file)+i)==col(file)
    val = sum(file[label])
    exp_val = val/length(file[label])
    file[label] = file[label]/exp_val
    file[lower.tri(file)] = t(file)[lower.tri(file)]
  }
  return(file)
}

con_data = obs_exp_cal(con_convert)
dif_data = obs_exp_cal(dif_convert)

# 画热图
pheatmap(con_data, cluster_rows=FALSE, cluster_cols=FALSE, border_color=NA,
         color=colorRampPalette(c("navy", "white", "firebrick3"))(50),
         filename="normal_zh.pdf")
pheatmap(dif_data, cluster_rows=FALSE, cluster_cols=FALSE, border_color=NA,
         color=colorRampPalette(c("navy", "white", "firebrick3"))(50),
         filename="abnormal_zh.pdf")
r1 = dif_data/con_data
pheatmap(r1, cluster_rows=FALSE, cluster_cols=FALSE, border_color=NA,
         color=colorRampPalette(c("navy", "white", "firebrick3"))(50),
         filename="abnormal/normal.pdf")
