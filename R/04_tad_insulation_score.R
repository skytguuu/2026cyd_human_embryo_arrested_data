# 4. TAD分析 - R部分: insulation score与DI

library("ggplot2")

# 读取insulation score数据
tad = read.table("../pn2_TAD_40kb.txt", header=TRUE)
filelist <- list.files(pattern=".bedGraph")
datalist <- lapply(filelist, function(x) read.table(x, header=FALSE, stringsAsFactors=F, skip=1))
m = length(filelist)

for (k in 1:m){
  data = as.data.frame(datalist[k])
  data[is.na(data)] = 0
  colnames(data) = c("Chr", "Start", "End", "Value")
  f_loc = c()
  f_value = c()
  for (i in 1:nrow(tad)){
    length = tad[i, 3] - tad[i, 2]
    up = 40000 * floor((tad[i, 2] - 0.5*length) / 40000)
    down = 40000 * ceiling((tad[i, 3] + 0.5*length) / 40000)
    part = subset(data, data$Chr == as.character(tad[i, 1]) & up <= data$Start & down > data$End)
    fgsc = c()
    fgsc_loc = c()
    for (j in 1:nrow(part)){
      if (part[j, 4] > 0.25){
        fgsc_loc = c(fgsc_loc, j*40000/length)
        fgsc = c(fgsc, part[j, 4])
      }
    }
    f_loc = c(f_loc, fgsc_loc)
    f_value = c(f_value, fgsc)
  }
  name = strsplit(filelist[k], split="\\.")
  cell = rep(name[[1]][1], length(f_loc))
  if (k == 1){
    mat = data.frame(location=f_loc, type=cell, value=f_value)
  } else {
    tmp = data.frame(location=f_loc, type=cell, value=f_value)
    mat = rbind(mat, tmp)
  }
}

r1 = ggplot(mat, aes(x=location, y=value)) + theme_bw() + theme(panel.grid=element_blank())
r2 = r1 + geom_smooth(data=mat[mat$type=="hum_8cell_2PN_40kb_input", ],
                       aes(color=type), size=2, method="loess",
                       color="black", span=0.05, se=FALSE, show.legend=TRUE)
r3 = r2 + geom_smooth(data=mat[mat$type=="hum_8cell_abnorm_2PN_40kb_input", ],
                       aes(color=type), size=2, method="loess",
                       span=0.05, se=FALSE)
r3

# DI
a1 = read.table("pn2_40kb_DI.txt", header=FALSE)
a2 = read.table("pn2_abnormal_40kb_DI.txt", header=FALSE)
boxplot(abs(a1$V4), abs(a2$V4), outline=FALSE, ylab="absolute DI",
        names=c("Normal", "Abnormal"))
