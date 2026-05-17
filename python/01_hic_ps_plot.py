# 1. 基本HiC数据分析 - Python部分: P(s)图

import warnings
warnings.filterwarnings("ignore")
from itertools import combinations
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
import pandas as pd
import bioframe
import cooler
import cooltools

# 设置染色体参数 - 人类
hg19_chromsizes = bioframe.fetch_chromsizes('hg19')
hg19_cens = bioframe.fetch_centromeres('hg19')
hg19_arms = bioframe.make_chromarms(hg19_chromsizes, hg19_cens)

# 读取文件并计算
con = cooler.Cooler('8cell_2PN_HiC.cool')
dif = cooler.Cooler('8cell_2PN_abnor_HiC.cool')
resolution = 500000

# 计算并画图
m = []
m1 = []
for i in [con, dif]:
    cvd = cooltools.expected_cis(clr=i, view_df=hg19_arms, smooth=True,
                                  aggregate_smoothed=True, nproc=8)
    cvd['s_bp'] = cvd['dist'] * resolution
    cvd['balanced.avg.smoothed'].loc[cvd['dist'] < 2] = np.nan
    cvd = cvd[cvd['s_bp'] >= 2*resolution]
    cvd_merged = cvd.drop_duplicates(subset=['dist'])[['s_bp', 'balanced.avg.smoothed.agg']]
    m.append([cvd_merged['s_bp'], cvd_merged['balanced.avg.smoothed.agg']])

    # 画着丝粒两端跨chromsome_arms的概率ps图
    cvd_inter = cooltools.expected_cis(clr=i, view_df=hg19_arms, intra_only=False)
    cvd_inter['s_bp'] = cvd_inter['dist'] * resolution
    m1.append([cvd_inter['s_bp'], cvd_merged['balanced.avg.smoothed.agg']])

# figure
f, ax = plt.subplots(1, 1)
p1, = ax.loglog(m[0][0], m[0][1], markersize=5, color="red", label="NC")
p2, = ax.loglog(m[1][0], m[1][1], markersize=5, color="orange", label="SH")
ax.set(xlabel='separation, bp', ylabel='IC contact frequency')
ax.grid(lw=0.5)
x = m[0][0]
r1, = plt.loglog(x, 1/(x*1e-6), color="black", linewidth=1, linestyle="-", label="S^(-1)")
r2, = plt.loglog(x, (x*1e-3)**-0.5, color="black", linewidth=1, linestyle="--", label="S^(-0.5)")
plt.legend()
# plt.savefig('hum_8cell_genome_contact_500kb.pdf')
plt.show()

# 画slope图
f, ax = plt.subplots(1, 1)
der1 = np.gradient(np.log(m[0][1]), np.log(m[0][0]))
der2 = np.gradient(np.log(m[1][1]), np.log(m[1][0]))
p1, = ax.semilogx(m[0][0], der1, alpha=0.5, markersize=5, color="red", label="Normal")
p2, = ax.semilogx(m[1][0], der2, alpha=0.5, markersize=5, color="orange", label="Abnormal")
ax.set(xlabel='P(s) derivative plot', ylabel='Slope')
ax.grid(lw=0.5)
plt.legend()
plt.savefig('hum8cell_ps_slope_500kb.pdf')
plt.show()
