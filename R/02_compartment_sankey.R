# 2. AB区室 - R部分: 桑吉图与AB转化比例

library("ggplot2")
library("ggalluvial")
library(dplyr)
library(reshape2)

d = read.table("cyd_chr1-22_mean_500kb.txt", header=TRUE, stringsAsFactors=FALSE)
d[is.na(d)] = 0

# A/B compartment分类
input = d[, 4:5]
input$NC = ifelse(d$hum_8cell_2PN_merge_500000_iced.matrix > 0, "A",
             ifelse(d$hum_8cell_2PN_merge_500000_iced.matrix < 0, "B", "Undefined"))
input$SH = ifelse(d$hum_8cell_abnorm_2PN_merge_500000_iced.matrix > 0, "A",
             ifelse(d$hum_8cell_abnorm_2PN_merge_500000_iced.matrix < 0, "B", "Undefined"))

# 转化为桑吉图格式
res = input %>% group_by(NC, SH) %>% summarise(count=n())
res$seq = rownames(res)
res1 = melt(as.data.frame(res), id.vars=c("seq", "count"),
            variable.name="cell", value.name="state")
res1$state = factor(res1$state, levels=c("A", "B", "Undefined"))
res1$cell = factor(res1$cell, levels=c("NC", "SH"))

p = ggplot(res1, aes(x=cell, stratum=state, alluvium=seq, y=count, fill=state, label=state))
p = p + scale_fill_manual(values=c("darkblue", "orange", "grey"))
p = p + geom_flow(aes(fill=state), stat="alluvium", lode.guidance="frontback", color="darkgray")
p = p + geom_stratum(color="white", size=2) + theme_bw() + theme(panel.grid=element_blank())
p

# 计算AB转化的比例
d = read.table("cyd_chr1-22_mean_500kb.txt", header=TRUE, stringsAsFactors=FALSE)
d[is.na(d)] = 0
nc_sh_a2b = subset(d, d$hum_8cell_2PN_merge_500000_iced.matrix > 0 &
                       d$hum_8cell_abnorm_2PN_merge_500000_iced.matrix < 0)
nc_sh_b2a = subset(d, d$hum_8cell_2PN_merge_500000_iced.matrix < 0 &
                       d$hum_8cell_abnorm_2PN_merge_500000_iced.matrix > 0)
nc_sh_stable = subset(d, d$hum_8cell_2PN_merge_500000_iced.matrix *
                          d$hum_8cell_abnorm_2PN_merge_500000_iced.matrix > 0)

data = data.frame(type=rep("normal2ab", 3),
                  group=c("A2B", "B2A", "Stable"),
                  count=c(dim(nc_sh_a2b)[1], dim(nc_sh_b2a)[1], dim(nc_sh_stable)[1]))

# 画柱状图
library(dplyr)
data %>%
  group_by(type) %>%
  mutate(prop=count/sum(count)) %>%
  ggplot(aes(x=type, y=count, fill=group)) +
  geom_bar(position="fill", stat="identity") +
  geom_text(aes(label=sprintf("%.1f%%", prop*100)),
            position=position_fill(vjust=0.5), size=3) +
  ggsci::scale_fill_nejm() +
  scale_y_continuous(labels=scales::percent) + theme_bw() +
  theme(legend.position="right")
