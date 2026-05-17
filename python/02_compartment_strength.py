# 2. AB区室 - Python部分: Compartment strength箱线图与统计

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon
import pandas as pd

# 画箱型图
nc = np.load('h8cell_normal_500kb_compartment_strength.saddledump.npz')
sh = np.load('h8cell_abnor_500kb_compartment_strength.saddledump.npz')
r1 = nc['saddle_strength']
r2 = sh['saddle_strength']

plt.figure(figsize=(10, 5))
plt.title('Compartment Strength', fontsize=20)
labels = 'Normal', 'Abnormal'
m = plt.boxplot([r1, r2], labels=labels, showfliers=False, patch_artist=True)
colors = ["#E69F00", "#56B4E9"]
for patch, color in zip(m['boxes'], colors):
    patch.set_facecolor(color)
# plt.savefig("compartment_strength_boxplot.pdf")
plt.show()

# python进行统计
stat, p = wilcoxon(r1, r2)

# 或者使用pandas保存为dataframe格式文件，用R来画图
df = pd.DataFrame([r1, r2], index=["NC", "SH"])
df.to_csv("compartment_saddle_strength.txt", sep="\t")
