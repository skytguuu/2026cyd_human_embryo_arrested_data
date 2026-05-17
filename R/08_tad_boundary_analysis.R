# 8. TAD boundary基因分析

library(GenomicRanges)
library(tidyverse)
library(ggpubr)

# 加载边界数据
norm = read.table("../tad/pn2_40kb_TAD_boundary.txt", header=FALSE)
abnor = read.table("../tad/pn2_abnromal_40kb_TAD_boundary.txt", header=FALSE)
norm_gr <- GRanges(norm$V1, IRanges(norm$V2, norm$V3))
abno_gr <- GRanges(abnor$V1, IRanges(abnor$V2, abnor$V3))

# 分类边界
conserved_idx <- findOverlaps(norm_gr, abno_gr)
conserved_gr <- norm_gr[unique(queryHits(conserved_idx))]
lost_gr <- norm_gr[-unique(queryHits(conserved_idx))]
gained_gr <- abno_gr[-unique(subjectHits(conserved_idx))]

boundary_list <- list(Conserved=conserved_gr,
                      Variable_Lost=lost_gr,
                      Variable_Gained=gained_gr)

# 准备边界分类数据
conserved_df <- as.data.frame(conserved_gr) %>% mutate(Type="Conserved")
lost_df <- as.data.frame(lost_gr) %>% mutate(Type="Variable (Lost)")
boundary_combined_df <- bind_rows(conserved_df, lost_df)
boundary_combined_gr <- GRanges(boundary_combined_df)

# 计算ZGA基因密度
buffer_size <- 20000
boundary_extended <- resize(boundary_combined_gr,
                            width=width(boundary_combined_gr) + 2*buffer_size,
                            fix="center")
boundary_combined_df$ZGA_Density <- countOverlaps(boundary_extended, zga_gr)

# 合并所有边界
df_cons <- as.data.frame(conserved_gr) %>% mutate(Group="Conserved")
df_lost <- as.data.frame(lost_gr) %>% mutate(Group="Lost")
df_gain <- as.data.frame(gained_gr) %>% mutate(Group="Gained")
boundary_all_df <- bind_rows(df_cons, df_lost, df_gain) %>%
  mutate(Group=factor(Group, levels=c("Conserved", "Lost", "Gained")))

all_gr <- GRanges(boundary_all_df)
all_gr_ext <- resize(all_gr, width=width(all_gr) + 80000, fix="center")
boundary_all_df$ZGA_Density <- countOverlaps(all_gr_ext, zga_gr)

my_cols <- c("Conserved"="#B3B3B3", "Lost"="#E64B35", "Gained"="#4DBBD5")
d1 = boundary_all_df[boundary_all_df$ZGA_Density != 0, ]

p_comparison <- ggplot(d1, aes(x=Group, y=ZGA_Density, fill=Group)) +
  geom_boxplot(outlier.shape=16, outlier.size=0.5, outlier.alpha=0.2, width=0.6) +
  scale_fill_manual(values=my_cols) +
  stat_compare_means(comparisons=list(c("Conserved", "Lost"), c("Conserved", "Gained")),
                     method="wilcox.test", label="p.signif") +
  theme_bw() +
  theme(panel.grid=element_blank(), legend.position="none",
        axis.text=element_text(size=12, color="black"),
        axis.title=element_text(size=13, face="bold")) +
  labs(x="TAD Boundary Dynamics (Arrested vs. Normal)",
       y="ZGA Gene Density (count per boundary)",
       title="Genomic Features of Variable TAD Boundaries")
p_comparison

# ZGA在可变边界的基因表达情况
a1 = read.csv("../rnaseq/hum_8cell_rnaseq_fpkm_symbol.csv", header=TRUE)
zga_gr1 <- GRanges(seqnames=zga1$chr,
                   ranges=IRanges(start=as.numeric(zga1$start), end=as.numeric(zga1$end)),
                   genes=zga1$symbol)
hits1 <- findOverlaps(gained_gr, zga_gr1)
hits2 <- findOverlaps(lost_gr, zga_gr1)
hits3 <- findOverlaps(conserved_gr, zga_gr1)
o1 <- zga_gr1$genes[subjectHits(hits1)]
o2 <- zga_gr1$genes[subjectHits(hits2)]
o3 <- zga_gr1$genes[subjectHits(hits3)]
da1 = a1[a1$symbol %in% o1, ]
da2 = a1[a1$symbol %in% o2, ]
da3 = a1[a1$symbol %in% o3, ]

# TE分析
te = read.table("hg19_TE", header=FALSE)
te_gr = GRanges(seqnames=te$V1, ranges=IRanges(start=te$V2, end=te$V3), gene=te$V4)

get_te_coverage <- function(gr, te_gr) {
  hits <- countOverlaps(gr, te_gr)
  return(hits / width(gr)[1] * 1000)
}

te_results <- data.frame(
  Group=rep(names(boundary_list), sapply(boundary_list, length)),
  TE_Density=unlist(lapply(boundary_list, get_te_coverage, te_gr=te_gr))
)

ggplot(te_results, aes(x=Group, y=TE_Density, fill=Group)) +
  geom_boxplot(outlier.shape=NA) +
  stat_compare_means(method="wilcox.test") + theme_bw() +
  labs(title="Transposable Element Enrichment at TAD Boundaries")

ggplot(te_results[te_results$Group != "Conserved", ],
       aes(x=factor(Group, levels=c("Variable_Lost", "Variable_Gained")),
           y=TE_Density, fill=Group)) +
  geom_boxplot(outlier.shape=NA) +
  stat_compare_means(comparisons=list(c("Variable_Lost", "Variable_Gained")),
                     method="wilcox.test", label="p.signif") + theme_bw() +
  labs(title="Transposable Element Enrichment at TAD Boundaries")
