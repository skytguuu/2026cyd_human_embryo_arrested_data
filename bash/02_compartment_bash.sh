#!/bin/bash
# 2. AB区室 - Bash部分

# 用cool来计算compartment
source activate cool5

# 计算AB区室 - 计算expected count
cooltools expected-cis 8cell_2PN_HiC.cool -o normal_obs_exp
cooltools eigs-cis 8cell_2PN_HiC.cool -o normal_cool5

# 计算compartment strength(signal), 公式是(AA+BB)/(AB+BA)
# 计算saddle并画图, strength文件保存在.npz文件中
cooltools saddle 8cell_2PN_HiC.cool normal_cool5.cis.vecs.tsv normal_obs_exp \
    --strength --qrange 0.2 0.8 \
    -o h8cell_normal_500kb_compartment_strength --fig pdf

# 进行这一步的原因，在于R输出的input格式尽管使用\t分隔，但是其中仍有空格的存在
# bedtools对input的要求为严格的\t分隔，所以没有这步会报错
ls *.txt | while read id; do sed -i 's/ //g' $id; done
ls *.txt | while read id; do
    perl ~/data1/script/perl_script/AB_compartment_merge.pl $id ${id%_chr*}_compart
    cat ${id%_chr*}_compart_A.bed ${id%_chr*}_compart_B.bed | sort -V -k1,1 -k2n,2 > ${id%_chr*}_compartment.bed
done

# AB区室基因 - 方法2：计算AB转化的区域
bedtools intersect -a hum_8cell_2PN_compart_A.bed -b hum_8cell_abnorm_2PN_compart_B.bed > h8cell_a2b.bed
perl ~/data1/script/perl_script/finf_compartment_genes.pl h8cell_a2b.bed \
    ~/data1/Annotation/RNAseq_anno/Genecode/hum/gencode_v19_gene_pos.txt h8cell_a2b.genes
