# 7. SSIM基因分析

library(GenomicRanges)
library(dplyr)
library(org.Hs.eg.db)
library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ggpubr)

load("20260323_ssim.RData")

safe <- results %>% filter(!is.na(start) & !is.na(end)) %>% filter(start >= 0 & end >= start)
zga = read.table("../rnaseq/hum_8vs4_ZGA_fc2_p0.01.txt", header=TRUE)
fpkm = read.table("../rnaseq/hum_8cell_merge_RNAseq1.txt", header=TRUE)

fpkm$ensembl_id <- sub("\\..*", "", fpkm$tracking_id)
fpkm$symbol <- mapIds(org.Hs.eg.db, keys=fpkm$ensembl_id,
                       column="SYMBOL", keytype="ENSEMBL", multiVals="first")
zga1 = fpkm[fpkm$symbol %in% zga$Symbol, ]
zga1$chr = stringr::str_split(zga1$locus, pattern=":", simplify=TRUE)[, 1]
zga1$loc = stringr::str_split(zga1$locus, pattern=":", simplify=TRUE)[, 2]
zga1$start = stringr::str_split(zga1$loc, pattern="-", simplify=TRUE)[, 1]
zga1$end = stringr::str_split(zga1$loc, pattern="-", simplify=TRUE)[, 2]

# ZGA基因密度
res_gr <- GRanges(seqnames=safe$chr,
                  ranges=IRanges(start=safe$start, end=safe$end),
                  group=safe$group, z_ssim=safe$z_ssim)
zga_gr <- GRanges(seqnames=zga1$chr,
                  ranges=IRanges(start=as.numeric(zga1$start), end=as.numeric(zga1$end)))

hits <- countOverlaps(res_gr, zga_gr)
safe$zga_density <- hits

safe1 = safe[!is.na(safe$group), ]
ggplot(safe1, aes(x=group, y=zga_density, fill=group)) +
  geom_boxplot(outlier.shape=NA) +
  stat_compare_means(method="wilcox.test") +
  theme_bw() +
  labs(title="ZGA Gene Density Enrichment", y="Number of ZGA Genes per Window")

# 比较AB区室强度
pc1 = read.table("../compartment/500kb/cyd_chr1-22_mean_500kb.txt",
                 header=TRUE, stringsAsFactors=FALSE)
pc1_gr <- GRanges(seqnames=pc1$seqnames,
                  ranges=IRanges(start=pc1$start, end=pc1$end),
                  normal=pc1$hum_8cell_2PN_merge_500000_iced.matrix,
                  abnor=pc1$hum_8cell_abnorm_2PN_merge_500000_iced.matrix)

hits <- findOverlaps(res_gr, pc1_gr)
overlap_df <- data.frame(query_idx=queryHits(hits),
                         normal_pc1=pc1_gr$normal[subjectHits(hits)],
                         abnor_pc1=pc1_gr$abnor[subjectHits(hits)])

intensity_res <- overlap_df %>%
  group_by(query_idx) %>%
  summarise(avg_normal_strength=mean(abs(normal_pc1), na.rm=TRUE),
            avg_abnor_strength=mean(abs(abnor_pc1), na.rm=TRUE))

results_final <- as.data.frame(res_gr) %>%
  mutate(query_idx=row_number()) %>%
  left_join(intensity_res, by="query_idx")

# 分面箱线图
plot_intensity <- results_final %>%
  filter(!is.na(avg_normal_strength)) %>%
  pivot_longer(cols=c(avg_normal_strength, avg_abnor_strength),
               names_to="Condition", values_to="Strength")
plot_intensity = plot_intensity[!is.na(plot_intensity$group.V1), ]

ggplot(plot_intensity, aes(x=group.V1, y=Strength, fill=group.V1)) +
  geom_boxplot(outlier.shape=NA, width=0.6) +
  facet_wrap(~Condition) +
  stat_compare_means(method="wilcox.test") +
  scale_fill_manual(values=c("sig"="#E64B35", "nodif"="#4DBBD5")) +
  theme_bw() +
  labs(title="Compartment Strength (abs(PC1))", y="Average Intensity",
       x="Structural Similarity")

# 合并Normal/Abnormal比较
plot_dat_combined <- results_final %>%
  filter(!is.na(avg_normal_strength) & !is.na(avg_abnor_strength)) %>%
  pivot_longer(cols=c(avg_normal_strength, avg_abnor_strength),
               names_to="Status", values_to="Strength") %>%
  mutate(Status=case_when(Status=="avg_normal_strength" ~ "Normal",
                           Status=="avg_abnor_strength" ~ "Abnormal"),
         Status=factor(Status, levels=c("Normal", "Abnormal")),
         group.V1=factor(group.V1, levels=c("nodif", "sig")))
plot_dat_combined = plot_dat_combined[!is.na(plot_dat_combined$group.V1), ]

ggplot(plot_dat_combined, aes(x=group.V1, y=Strength, fill=Status)) +
  geom_boxplot(outlier.shape=NA, width=0.7, position=position_dodge(0.8)) +
  scale_fill_manual(values=c("Normal"="#74add1", "Abnormal"="#f46d43")) +
  stat_compare_means(aes(group=Status), label="p.signif", method="wilcox.test") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text=element_text(color="black", size=10),
        legend.position="top") +
  labs(x="Structural Similarity Group", fill="Condition")
