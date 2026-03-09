import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import os
from datetime import datetime

# ==========================================
# 1. 读取原始完整数据
# ==========================================
# 请确保 'all_reads_compare.csv' 文件与此脚本在同一目录下
df = pd.read_csv('all_reads_compare.csv')

# 处理列名可能存在的问题，确保前三列分别是 Gene, NC, dCas9
if 'Unnamed: 0' in df.columns:
    df.rename(columns={'Unnamed: 0': 'Gene'}, inplace=True)
elif df.columns[0] == '':
    df.columns.values[0] = 'Gene'
elif 'Gene' not in df.columns:
    df.columns = ['Gene', 'NC', 'dCas9']

df = df.dropna()

# ==========================================
# 2. 数据修剪 (Data Trimming) - 关键步骤
# ==========================================
# 为了在线性轴上清晰显示大多数数据，去除掉 dCas9 组极少数高得离谱的“超级富集”点
# 这里设定阈值为 dCas9 的第 95 百分位数 (即去除 top 5% 的极端数据)
upper_limit = df['dCas9'].quantile(0.95)
df_trimmed = df[df['dCas9'] <= upper_limit].copy()

# ==========================================
# 3. 数据转换与准备 (基于修剪后的数据，保持线性轴)
# ==========================================
# --- Panel A (小提琴图) 的数据准备 ---
df_melt = df_trimmed.melt(id_vars='Gene', value_vars=['NC', 'dCas9'], 
                          var_name='Group', value_name='Reads')
# 仅修改图中显示名称，不影响原始数据列名
df_melt['Group'] = df_melt['Group'].replace({'dCas9': 'WORF'})

# --- Panel B (热图) 的数据准备 ---
# 计算 dCas9 和 NC 之间的绝对 Reads 差异，用于排序
df_trimmed['Diff'] = df_trimmed['dCas9'] - df_trimmed['NC']

# 按照绝对差异降序排列，让富集效果最好的区域排在最上面
df_sorted = df_trimmed.sort_values(by='Diff', ascending=False).reset_index(drop=True)

# 准备热图专属数据矩阵 (只提取原始 Reads 列)
heat_df = df_sorted[['Gene', 'NC', 'dCas9']].set_index('Gene')
heat_df.columns = ['NC', 'WORF']

# ==========================================
# 4. 全局排版参数设置 (面上申请书/期刊通用标准)
# ==========================================
mpl.rcParams['font.family'] = 'Arial'        # 采用无衬线字体Arial
mpl.rcParams['pdf.fonttype'] = 42            # 确保保存的PDF在AI中字体可编辑
mpl.rcParams['axes.linewidth'] = 1.0         # 轴线宽度

# 创建一个 1行2列 的画板，Panel A 与 Panel B 宽度比为 2:1
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 5), dpi=300, 
                               gridspec_kw={'width_ratios': [2, 1], 'wspace': 0.35})

# ==========================================
# 5. 绘制 Panel A (小提琴图 - 线性轴)
# ==========================================
# 自定义颜色 (浅蓝对浅红)
colors = {'NC': '#A0CBE8', 'WORF': '#FF9D9A'}

# 画小提琴主体
sns.violinplot(x='Group', y='Reads', data=df_melt, 
               palette=colors, inner=None, linewidth=0, alpha=0.6, ax=ax1, zorder=1)

# 叠加箱线图 (统计学分布)
sns.boxplot(x='Group', y='Reads', data=df_melt, 
            width=0.15,
            boxprops={'facecolor':'none', 'edgecolor':'black', 'zorder':2, 'linewidth':0.5},
            medianprops={'color':'black', 'linewidth':1, 'zorder':2},
            whiskerprops={'linewidth':0.5, 'color':'black'}, 
            capprops={'linewidth':0.5, 'color':'black'},
            showfliers=False, ax=ax1)

# 叠加散点 (真实靶点的位置)
sns.stripplot(x='Group', y='Reads', data=df_melt, 
              palette=colors, size=2.5, jitter=0.15, alpha=0.8, ax=ax1, zorder=3, edgecolor='none')

# Panel A 的细节格式化
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.set_xlabel('')
ax1.set_ylabel('Reads', fontsize=12) # 原始Reads标签
ax1.tick_params(axis='both', which='major', labelsize=11)

# ==========================================
# 6. 绘制 Panel B (热图 - 线性颜色映射)
# ==========================================
# 使用渐变色谱 (从浅黄过渡到深红近黑)
cmap = sns.color_palette("rocket_r", as_cmap=True)

# 绘制热图
sns.heatmap(heat_df, 
            cmap=cmap, 
            ax=ax2, 
            yticklabels=False,  # 隐藏具体基因名，突出富集趋势
            cbar_kws={'label': 'Reads', 'shrink': 0.7, 'aspect': 30}, 
            linewidths=0)       # 去除白色分割线

# Panel B 的细节格式化
ax2.set_ylabel('Target Regions (Ranked by Absolute Enrichment)', fontsize=12)
ax2.tick_params(axis='x', which='major', labelsize=11, rotation=0)

# 调整 Colorbar 字体大小
cbar = ax2.collections[0].colorbar
cbar.ax.tick_params(labelsize=9)
cbar.set_label('Reads', size=11)

# ==========================================
# 7. 保存与输出
# ==========================================
# 保存为可编辑矢量图SVG，输出目录为 /output，文件名为 Fig6+时间戳
output_dir = os.path.join(os.getcwd(), 'output')
os.makedirs(output_dir, exist_ok=True)
timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
output_file = os.path.join(output_dir, f'Fig6+{timestamp}.svg')
plt.savefig(output_file, bbox_inches='tight', format='svg')

print(f"已剔除Top 5%的极端异常值，图表已保存为 {output_file}")