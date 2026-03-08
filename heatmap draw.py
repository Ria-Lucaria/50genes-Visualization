import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

# 1. 读取原始完整数据
# 请确保 'all_reads_compare.csv' 文件与此脚本在同一目录下
df = pd.read_csv('all_reads_compare.csv')
df.columns = ['Gene', 'NC', 'dCas9']
df = df.dropna()

# 2. 数据转换与计算
# 计算 log2(Reads + 1) 以便在热图中展示，避免被极端高值掩盖低表达基因
df['log2_NC'] = np.log2(df['NC'] + 1)
df['log2_dCas9'] = np.log2(df['dCas9'] + 1)

# 计算上调倍数 Log2 Fold Change (dCas9相对于NC的倍数变化)
df['Log2FC'] = df['log2_dCas9'] - df['log2_NC']

# 3. 数据排序
# 按照 Log2FC 从大到小降序排列，打造视觉上平滑的富集梯度
df_sorted = df.sort_values(by='Log2FC', ascending=False).reset_index(drop=True)

# 准备专门用于绘制热图的数据矩阵
heat_df = df_sorted[['Gene', 'log2_NC', 'log2_dCas9']].set_index('Gene')
heat_df.columns = ['NC', 'dCas9']

# ==========================================
# 4. 《Nature》期刊标准排版参数设置 (全局)
# ==========================================
mpl.rcParams['font.family'] = 'Arial'          # 强制使用 Arial 无衬线字体
mpl.rcParams['font.size'] = 7                  # 全局基础字号 7pt
mpl.rcParams['axes.linewidth'] = 0.5           # 坐标轴线宽 0.5pt
mpl.rcParams['pdf.fonttype'] = 42              # 保证导出 PDF 时字体为 TrueType (可编辑)

# 5. 创建画布
# 宽度 1.8 英寸 (约 45 mm，适合单栏拼图的侧边栏)，高度 6.0 英寸以容纳 52 个基因名
fig, ax = plt.subplots(figsize=(1.8, 6.0), dpi=300)

# 选择色盲友好且感知均匀的颜色渐变板 (从浅黄过渡到深红再到近黑)
cmap = sns.color_palette("rocket_r", as_cmap=True)

# 6. 绘制热图
sns.heatmap(heat_df, 
            cmap=cmap, 
            ax=ax, 
            yticklabels=True,                                        # 强制显示所有 Y 轴标签 (基因名)
            cbar_kws={'label': '$\log_2$(Reads + 1)', 'shrink': 0.6, 'aspect': 30}, 
            linewidths=0.2,                                          # 色块之间加入 0.2pt 的纯白分割线
            linecolor='white')

# 7. 细节打磨与格式化
# Y 轴 (基因名称) 格式化
ax.set_ylabel('')                                      # 隐藏多余的 "Gene" 标题
ax.tick_params(axis='y', which='major', labelsize=6, rotation=0) # 基因名设置为 6pt，不旋转
ax.yaxis.set_tick_params(pad=2)                        # 微调基因名与色块的间距

# X 轴 (NC / dCas9 标签) 格式化
ax.tick_params(axis='x', which='major', labelsize=7, rotation=45) # 组别名称 7pt，旋转45度防重叠

# 图例 (Colorbar) 格式化
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=5, length=2, width=0.5)  # Colorbar 刻度数字 5pt，刻度线 0.5pt
cbar.set_label('$\log_2$(Reads + 1)', size=6)          # Colorbar 标题 6pt
cbar.outline.set_linewidth(0.5)                        # Colorbar 外框线宽 0.5pt

# 8. 自动调整布局并导出矢量图与高清位图
plt.tight_layout()
plt.savefig('full_52genes_heatmap_nature.png', dpi=300, bbox_inches='tight')
plt.savefig('full_52genes_heatmap_nature.pdf', bbox_inches='tight') # 投递 Nature 等期刊首选该格式
plt.close()

print("Full 52-gene heatmap has been successfully saved!")