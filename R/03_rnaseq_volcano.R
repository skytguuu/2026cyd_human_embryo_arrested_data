# 3. RNAseq分析 - 火山图

library(ggplot2)
library(ggrepel)

res$change <- as.factor(
    ifelse(
        res$adj.P.Val < 0.05 & abs(res$logFC) > log2(4),
        ifelse(res$logFC > log2(4), "Up", "Down"),
        "NoDiff"
    )
)

valcano <- ggplot(data=res, aes(x=logFC, y=-log10(adj.P.Val), color=change)) +
  geom_point(alpha=0.8, size=1) +
  theme_bw(base_size=15) +
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  ggtitle("Valcano") +
  scale_color_manual(name="", values=c("red", "green", "black"),
                     limits=c("Up", "Down", "NoDiff")) +
  geom_vline(xintercept=c(-log2(4), log2(4)), lty=2, col="gray", lwd=0.5) +
  geom_hline(yintercept=-log10(0.05), lty=2, col="gray", lwd=0.5)
valcano

# 标记基因 - 上下调各10个
res$label = ''
res <- res[order(res$adj.P.Val), ]
up <- head(rownames(res)[which(res$change=="Up")], 10)
down <- head(rownames(res)[which(res$change=="Down")], 10)
top10 <- c(as.character(up), as.character(down))
res$label[match(top10, rownames(res))] <- top10

# 画图
p <- ggplot(res, aes(x=logFC, y=-log10(adj.P.Val))) +
  geom_hline(aes(yintercept=-log10(0.05)), color="#999999", linetype="dashed", linewidth=1) +
  geom_vline(aes(xintercept=-2), color="#999999", linetype="dashed", linewidth=1) +
  geom_vline(aes(xintercept=2), color="#999999", linetype="dashed", linewidth=1) +
  geom_point(color="grey", alpha=0.5) +
  geom_point(data=res[res$adj.P.Val < 0.05 & (res$logFC > 2), ],
             stroke=0.5, size=2, shape=16, color="firebrick", alpha=0.5) +
  geom_point(data=res[res$adj.P.Val < 0.05 & (res$logFC < -2), ],
             stroke=0.5, size=2, shape=16, color="dodgerblue", alpha=0.5) +
  labs(x="Log2 fold change", y="-Log10(pvalue)") +
  theme_classic() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.title=element_text(size=12),
        axis.text=element_text(size=12, color="black"),
        plot.title=element_text(hjust=0.5, size=12)) +
  theme(legend.position="none") +
  geom_text_repel(aes(label=`label`), color="black", size=5, fontface="italic",
                  arrow=arrow(ends="first", length=unit(0.01, "npc")),
                  box.padding=0.2, point.padding=0.3,
                  segment.color='black', segment.size=0.3, force=1,
                  max.iter=3e3, max.overlaps=300) +
  xlim(-100, 100) + ylim(0, 10)
