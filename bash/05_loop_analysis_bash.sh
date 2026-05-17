#!/bin/bash
# 5. Loop分析 - Bash部分

# 统计每个染色体的染色质环个数差异
# 因为chr16 -520000 -480000有负值，会报错，去掉
awk '$2 >= 0 && $3 > $2' pn2_enh_1e-4.bed > output1.txt
bedtools intersect -a output1.txt -b output2.txt -wa | uniq | wc -l
bedtools intersect -a output1.txt -b output2.txt -wa | uniq > nor_uniq_loop.txt

# loop画图
source activate cool5
hicAggregateContacts --matrix ../tad/40kb/hum_8cell_2PN_merge_40000.cool \
    --BED output1.txt --outFileName normal_loop_Contacts \
    --range 400000:2000000 --numberOfBins 10 --operationType mean \
    --transform obs/exp --outFileContactPairs normal_Contacts --mode intra-chr

# 转录因子分析
findMotifsGenome.pl nor_uniq_loop.txt hg19 nor_uniq_motif/ -len 8,10,12

# 差异基因和染色质环的转录因子
grep -f nor_spe_deg.txt nor_uniq_loop.txt > 1.txt
grep -f abnor_spe_deg.txt abnor_uniq_loop.txt > 2.txt
findMotifsGenome.pl 1.txt hg19 nor_uniq_deg_motif/ -len 8,10,12
