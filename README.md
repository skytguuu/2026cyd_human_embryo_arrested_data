# Hi-C数据分析代码

曹玥迪Hi-C数据分析项目的代码集合，包含Bash、R和Python三种语言的脚本。

## 目录结构

```
HiC_analysis/
├── bash/                          # Bash脚本
│   ├── 01_hic_basic_analysis.sh       # 基本HiC数据分析（.hic文件生成、cis/trans比例、cool文件转换）
│   ├── 02_compartment_bash.sh         # AB区室计算（cooltools计算compartment、saddle）
│   ├── 04_tad_analysis_bash.sh        # TAD分析（insulation score、average TAD、RTI、boundary）
│   ├── 05_loop_analysis_bash.sh       # Loop分析（loop统计、转录因子motif分析）
│   └── 06_sv_analysis.sh             # SV分析（CNV、EagleC预测SV、hic_breakfinder、HiNT）
│
├── R/                             # R脚本
│   ├── 01_hic_basic_cis_trans.R       # cis/trans比例箱线图
│   ├── 01_hic_inter_chromosome_heatmap.R  # 跨染色体互作热图
│   ├── 01_hic_inter_chr_scatter.R    # 跨染色体散点图与比例柱状图
│   ├── 02_compartment_boxplot.R      # Compartment strength箱线图
│   ├── 02_compartment_pca.R          # 使用HiTC计算AB区室
│   ├── 02_compartment_sankey.R       # 桑吉图与AB转化比例
│   ├── 02_compartment_density.R      # PC密度图与散点密度图
│   ├── 02_compartment_gene_expression.R  # AB区室基因与表达分析
│   ├── 02_compartment_a2b_go.R       # AB转化区域GO富集分析
│   ├── 03_rnaseq_deg.R              # RNAseq差异基因与通路富集
│   ├── 03_rnaseq_updown_go.R        # 上下调基因GO与GSEA
│   ├── 03_rnaseq_volcano.R          # 火山图
│   ├── 04_tad_insulation_score.R     # TAD insulation score与DI
│   ├── 04_tad_rti_boundary.R        # TAD RTI与边界分析
│   ├── 04_tad_enrichment.R          # TAD富集分析
│   ├── 05_loop_analysis.R           # Loop分析（GO富集、柱状图）
│   ├── 05_loop_ab_compartment.R     # AB区室转化的loop基因
│   ├── 06_mdecay_analysis.R         # M-decay基因分析
│   ├── 06_minor_zga.R              # Minor ZGA分析
│   ├── 07_ssim_analysis.R          # SSIM基因分析
│   └── 08_tad_boundary_analysis.R   # TAD boundary基因分析
│
└── python/                        # Python脚本
    ├── 01_hic_ps_plot.py             # P(s)图与slope图
    └── 02_compartment_strength.py    # Compartment strength箱线图与统计
```

## 分析内容

1. **基本HiC数据分析** - .hic文件生成、cis/trans比例、P(s)图、跨染色体互作
2. **AB区室分析** - compartment计算、saddle图、桑吉图、AB转化基因与GO富集
3. **RNAseq分析** - 差异基因(limma)、GO/GSEA富集、火山图
4. **TAD分析** - insulation score、average TAD、RTI、boundary、TAD富集
5. **Loop分析** - 染色质环统计、转录因子分析、AB转化loop基因
6. **SV分析** - CNV、EagleC预测SV、hic_breakfinder、HiNT
7. **M-decay基因分析** - Z-decay/M-decay分类、minor ZGA
8. **SSIM基因分析** - 结构相似性与ZGA密度、compartment强度
9. **TAD boundary基因分析** - 保守/可变边界分类、TE富集

## 依赖环境

- **Bash**: HiC-Pro, cooltools, bedtools, cooler, HiCExplorer, coolpup.py
- **R**: HiTC, ggplot2, ggpubr, ggalluvial, pheatmap, clusterProfiler, limma, GenomicRanges, org.Hs.eg.db
- **Python**: cooler, cooltools, bioframe, matplotlib, numpy, pandas, scipy
