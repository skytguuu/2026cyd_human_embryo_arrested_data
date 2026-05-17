#!/bin/bash
# 6. SV分析 - Bash部分

# 增强子劫持 - neoloopfinder工具
# 计算CNV
calculate-cnv -H ../tad/40kb/hum_8cell_2PN_merge_40000.cool -g hg19 -e DpnII \
    --output normal_40kb.CNV-profile.bedGraph
segment-cnv --cnv-file normal_40kb.CNV-profile.bedGraph --binsize 40000 --ploidy 2 \
    --output normal_40kb.CNV-seg.bedGraph --nproc 12
plot-cnv --cnv-profile normal_40kb.CNV-profile.bedGraph \
    --cnv-segment normal_40kb.CNV-seg.bedGraph \
    --output-figure-name normal_40k.CNV.genome-wide.png \
    --dot-size 0.5 --dot-alpha 0.2 --line-width 1 --boundary-width 0.5 \
    --label-size 7 --tick-label-size 6 --clean-mode

# 矫正Hi-C通过CNV
correct-cnv -H ../tad/40kb/hum_8cell_2PN_merge_40000.cool \
    --cnv-file normal_40k.CNV-seg.bedGraph --nproc 12 -f

# 用EagleC软件来预测SV
predictSV --hic-5k ../cool/normal_10kb.cool --hic-10k ../cool/normal_10kb.cool \
    --hic-50k ../cool/normal_50kb.cool -O normal -g hg19 --balance-type ICE \
    --output-format full --prob-cutoff-5k 0.8 --prob-cutoff-10k 0.8 --prob-cutoff-50k 0.99999

annotate-gene-fusion --sv-file abnormal.CNN_SVs.5K_combined.txt \
    --output-file abnormal.gene-fusions.txt \
    --buff-size 10000 --skip-rows 1 --ensembl-release 93 --species human

# 继续neoloopfinder
assemble-complexSVs -O abnormal -B abnormal.CNN_SVs.5K_combined.txt \
    --balance-type ICE --protocol insitu --nproc 24 \
    -H ../cool/abnormal_10kb.cool ../cool/abnormal_5kb.cool ../cool/abnormal_50kb.cool

# 安装hic_breakfinder
# 1. 先安装依赖eigen和bamtools
# 安装eigen
source activate hic
mamba install -c ccordoba12 eigen3

# 安装bamtools，不能使用conda安装bamtools
# 下载bamtools-2.5.1.zip
# 优先安装cmake
mamba install -c conda-forge/label/gcc7 cmake

# 在bamtools文件夹下
mkdir build
cd build
# 注意GCC版本，g++ --version，版本是11，所以需要设置为11
cmake -DCMAKE_CXX_STANDARD=11 -DCMAKE_INSTALL_PREFIX=/home/dell/data1/softwares/HiC/bamtools-2.5.1_master ..
make
make install

# 2. 安装hic_breakfinder
git clone https://github.com/dixonlab/hic_breakfinder.git
./configure --prefix=/home/dell/data1/softwares/HiC/hic_breakfinder \
    CPPFLAGS="-I /home/dell/data1/softwares/HiC/bamtools-2.5.1_master/include -I /home/dell/anaconda3/envs/hic/include/eigen3" \
    LDFLAGS="-L/home/dell/data1/softwares/HiC/bamtools-2.5.1_master/lib"
cmake CXXFLAGS += -I/home/dell/data1/softwares/HiC/bamtools-2.5.1/src/api
make
make install

# 超算使用hic_breakfinder
#SBATCH --job-name=hic_breakfinder
#SBATCH --partition=small
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH -N 1
#SBATCH --ntasks-per-node=5
module load hic_breakfinder/1.0-gcc-9.2.0
hic_breakfinder --bam-file ${Example}/K562_in_house_b38d5.nodup.bam \
    --exp-file-inter ${Example}/inter_expect_1Mb.hg38.txt \
    --exp-file-intra ${Example}/intra_expect_100kb.hg38.txt \
    --min-1kb --name Out

# 3. HINT软件找SV
source activate hinc
hint tl -m ../../cyd/hum_8cell_abnorm_2PN_merge.allValidPairs.hic \
    -f juicer --refdir ~/data1/Annotation/HiC_needfile/hic_breaker/hg19 \
    --backdir ~/data1/Annotation/HiC_needfile/hic_breaker/background_hg19/ \
    -g hg19 -n B1 -c 0.05 \
    --ppath /home/dell/anaconda3/envs/hinc/bin/pairix -p 12 \
    -o abnormal_HiNTtransl_juicerOUTPUT
