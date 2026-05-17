#!/bin/bash
# 4. TAD分析 - Bash部分

# 计算TAD的insulation score
ls *40000_iced.matrix | while read id; do
    python ~/data1/softwares/HiC/HiC-Pro_2.11.1/bin/utils/sparseToDense.py \
        --bins ~/data1/Annotation/HiC_needfile/human_genome_40kb.bed $id \
        -o ${id%_merge*}_40kb_sparse.matrix
done

ls *40kb_sparse.matrix | while read id; do
    paste ~/data1/Annotation/HiC_needfile/hum_40kb_ins_row.txt $id | \
        cat ~/data1/Annotation/HiC_needfile/hum_40kb_ins_col.txt - > ${id%_sparse*}_input.matrix
done

source activate hic
ls *input.matrix | while read id; do
    (nohup perl ~/data1/softwares/HiC/crane-nature-2015/scripts/matrix2insulation.pl \
        -i $id -b 2000000 -ids 400000 -im mean -bmoe 3 -nt 0.1 -v \
        1> $id.log &)
done

# 计算TAD的average TAD
ls *_500000.matrix | while read id; do
    (nohup hicConvertFormat -m $id \
        --bedFileHicpro ~/data1/Annotation/HiC_needfile/human_genome_500kb.bed \
        --inputFormat hicpro --outputFormat cool \
        -o ${id%.*}.cool >${id%.*}.log &)
done
# normalization
ls *.cool | while read id; do cooler balance $id; done

source activate cooltools
coolpup.py hum_8cell_2PN_merge_40000.cool ../pn2_TAD_40kb.txt \
    --features-format bed --rescale --local --outname normal_average_40TAD
plotpup.py --input normal_average_40TAD --scale log --output h8cell_normal.pdf

# TAD RTI
nohup Rscript ~/data1/script/Rscript/script/final_tad_RTI_score_args.r \
    ~/data1/Annotation/HiC_needfile/genome_40kb.bed SH_TAD_40kb.txt \
    ../../NC_SH_result/SH_merge_40000_iced.matrix SH_RTI_40kb_TAD --no-save

# TAD boundary - 定义小于400kb的为boundary
sort -V -k1,1 -k2n,2 NC_TAD_40kb.txt | \
    bedtools complement -i - -g ~/data1/Annotation/Genome_anno/chrom_mm9.sizes | \
    awk '{if(($3-$2)<=400000)print $1"\t"$2"\t"$3}' > NC_40kb_TAD_boundary.txt

# tad交集
bedtools intersect -a pn2_TAD_40kb.txt -b pn2_abnormal_TAD_40kb.txt -wa | uniq | wc -l
bedtools intersect -a pn2_abnormal_TAD_40kb.txt -b pn2_TAD_40kb.txt -wa | uniq | wc -l

# tad内部的基因
perl ~/data1/script/perl_script/allgene_located_compartment.pl pn2_TAD_40kb.txt \
    ~/data1/softwares/HiC/PSYCHIC-master/examples/hg19.genes.bed normal_TAD_40kb_genes.txt
perl ~/data1/script/perl_script/allgene_located_compartment.pl pn2_40kb_TAD_boundary.txt \
    ~/data1/softwares/HiC/PSYCHIC-master/examples/hg19.genes.bed normal_TAD_boundary_genes.txt
