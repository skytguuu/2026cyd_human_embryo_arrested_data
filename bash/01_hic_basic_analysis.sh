#!/bin/bash
# 1. 基本HiC数据分析 - Bash部分

# 生成.hic文件
nohup ~/data1/softwares/tools/Hi-C/HiC-Pro_2.11.1/bin/utils/hicpro2juicebox.sh \
    -i GSE135104_3.allValidPairs \
    -g ~/softwares/tools/Hi-C/HiC-Pro_2.11.1/annotation/chrom_mm9.sizes \
    -j ~/softwares/tools/Hi-C/juicer/juicer_tools.1.8.9_jcuda.0.8.jar \
    -r ~/softwares/tools/Hi-C/HiC-Pro_2.11.1/annotation/mm9_MboI.bed &

# 图40kb折线图 - 统计规定范围的cis/trans比例
ls *.allValidPairs | while read id; do
    for i in {40000,80000,100000,120000,2000000}; do
        perl ~/data1/script/perl_script/defined_dis_sum_ratio.pl $id $i
    done
done

# 计算<2Mb和>2Mb的interaction比例
# perl ~/data1/script/perl_script/cal_distance_hic_ratio.pl NC_merge.allValidPairs 2000000 19 > NC_2mb_ratio.txt
ls *.allValidPairs | while read id; do
    (nohup perl ~/data1/script/perl_script/cal_distance_hic_ratio.pl $id 2000000 19 > $id.txt &)
done
# 小鼠是19条，人是22条常染色体

# Ps图 - 格式转化，转化为cool文件
ls *_500000.matrix | while read id; do
    (nohup hicConvertFormat -m $id \
        --bedFileHicpro ~/data1/Annotation/HiC_needfile/human_genome_500kb.bed \
        --inputFormat hicpro --outputFormat cool \
        -o ${id%.*}.cool >${id%.*}.log &)
done

# normalization
ls *.cool | while read id; do cooler balance $id; done

# 跨染色体互作图 - obs/exp chr
perl ~/data1/script/perl_script/static_allvalid_chromsomeall_reads.pl ../NC/NC/NC_merge.allValidPairs | \
    grep -v "random" | grep -v "chrM" | grep -v "chrY" | grep -v "chrX" | \
    sort -V -k1,1 -k2,2 | sed '$d' > NC_zh.txt

perl ~/data1/script/perl_script/static_allvalid_chromsomeall_reads.pl ../SH/SH/SH_merge.allValidPairs | \
    grep -v "random" | grep -v "chrM" | grep -v "chrY" | grep -v "chrX" | \
    sort -V -k1,1 -k2,2 | sed '$d' > SH_zh.txt

# 统计单一染色体内部与其它染色体的cis/trans比例
ls *.allValidPairs | while read id; do
    perl ~/data1/script/perl_script/static_allvalid_chromsomeall_reads.pl $id | \
        sed '1d' | grep -v "random" | grep -v "chrM" | grep -v "chrY" | grep -v "chrX" | grep -v "chrUn" | \
        sort -V -k1,1 -k2,2 > ${id%.*}_inter_chr.txt
done

# 统计GC content和Gene density
# 1. bin基因组(1Mb)
cooltools genome binnify ~/data1/Annotation/Genome_anno/chrom_hg19.sizes 1000000 > hum_bins1mb.txt
# 2. 计算范围内的gc含量
cooltools genome gc hum_bins1mb.txt ~/data1/Annotation/fasta/hg19.fa > hg19_1mb_gc.txt
# 3. 计算范围内的gene density
cooltools genome genecov hum_bins1mb.txt hg19 > hg19_1mb_genomedensity.txt

# bedtools计算overlap
bedtools intersect -a ~/Annotation/Chip_anno/mm9_1mb_bin.bed -b h_ICM_TAD.txt -wa -wb > TAD_1mb.txt
# bedtools统计1mb范围内的个数
bedtools groupby -i zyb_control_TAD_1mb.txt -g 1,2,3 -c 4 -o count > zyb_control_TAD_1mb_num.txt
# 合并一句命令
bedtools intersect -a ~/Annotation/Chip_anno/mm9_1mb_bin.bed -b zyb_control_20kb_TAD.txt -wa -wb | \
    bedtools groupby -i - -g 1,2,3 -c 4 -o count > zyb_control_TAD_1mb.txt
