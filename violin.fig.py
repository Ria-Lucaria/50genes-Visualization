import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

# ==========================================
# 1. 读取原始完整数据
# ==========================================
# 请确保 'all_reads_compare.csv' 文件与此脚本在同一目录下
df = pd.read_csv('all_reads_compare.csv')
df.columns = ['Gene', 'NC', 'dCas9']
df = df.dropna()

# ==========================================
# 2. 数据转换 (Data Transformation)
# ==========================================
# 因为原始数据跨度极大（从个位数到上万），必须进行 Log10 转换，否则图形会被极端值压扁
df['log10_NC'] = np.log10(df['NC'] + 1)
df['log10_dCas9'] = np.log10(df['dCas9'] + 1)

# 将宽表 (Wide format) 转换为长表 (Long format)，以适应 Seaborn 的绘图逻辑
df_melt = df.melt(id_vars='Gene', value_vars=['log10_NC', 'log10_dCas9'], 
                  var_name='Group', value_name='log10_Reads')

# 清理标签名称，使其在图表 X 轴上显示得更简洁
df_melt['Group'] = df_melt['Group'].replace({'log10_NC': 'NC', 'log10_dCas9': 'dCas9'})

# ==========================================
# 3. 《Nature》期刊标准排版参数设置 (全局)
# ==========================================
mpl.rcParams['font.family'] = 'Arial'          # 强制使用 Arial 无衬线字体
mpl.rcParams['font.size'] = 7                  # 全局基础字号 7pt
mpl.rcParams['axes.linewidth'] = 0.5           # 坐标轴主线宽 0.5pt
mpl.rcParams['xtick.major.width'] = 0.5        # X轴刻度线宽 0.5pt
mpl.rcParams['ytick.major.width'] = 0.5        # Y轴刻度线宽 0.5pt
mpl.rcParams['xtick.direction'] = 'out'        # X轴刻度朝外
mpl.rcParams['ytick.direction'] = 'out'        # Y轴刻度朝外
mpl.rcParams['pdf.fonttype'] = 42              # 保证导出 PDF 时字体为 TrueType (可编辑)

# ==========================================
# 4. 核心绘图逻辑 (Plotting)
# ==========================================
# 创建画布：宽度 1.8 英寸，高度 2.5 英寸 (极其精巧的 Panel 尺寸)
fig, ax = plt.subplots(figsize=(1.8, 2.5), dpi=300)

# 设定色盲友好型对比色：NC (蓝色), dCas9 (朱红色)
colors = {'NC': '#56B4E9', 'dCas9': '#D55E00'}

# 层次1：绘制底层小提琴图 (概率密度轮廓)
sns.violinplot(x='Group', y='log10_Reads', data=df_melt, 
               palette=colors, inner=None, linewidth=0.5, ax=ax, zorder=1)

# 将小提琴的面色设置为 30% 透明度，以凸显内部的箱线和散点
for collection in ax.collections:
    collection.set_alpha(0.3)

# 层次2：叠加极简风格的内部箱线图 (展示中位数、四分位距)
sns.boxplot(x='Group', y='log10_Reads', data=df_melt, 
            width=0.15,                                              # 极窄的箱体
            boxprops={'facecolor':'none', 'edgecolor':'black', 'zorder':2, 'linewidth':0.5},
            medianprops={'color':'black', 'linewidth':1, 'zorder':2}, # 加粗黑色中位数线
            whiskerprops={'linewidth':0.5, 'color':'black'}, 
            capprops={'linewidth':0.5, 'color':'black'},
            showfliers=False, ax=ax)                                 # 不显示离群点（由散点图负责）

# 层次3：叠加随机抖动散点 (展示每一个真实靶点的位置)
sns.stripplot(x='Group', y='log10_Reads', data=df_melt, 
              palette=colors, size=2.5, jitter=0.15, alpha=0.8, ax=ax, zorder=3, edgecolor='none')

# ==========================================
# 5. 细节打磨与格式化
# ==========================================
# 去除右侧和顶部的边框 (Top & Right Spines)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# 设置轴标签
ax.set_xlabel('')                                           # 底部标签已明确是 NC/dCas9，无需 "Group" 字样
ax.set_ylabel('$\log_{10}(\mathrm{Reads} + 1)$', fontsize=7) # 标准数学符号标注

# 格式化刻度字体
ax.tick_params(axis='y', which='major', labelsize=6)
ax.tick_params(axis='x', which='major', labelsize=7)

# ==========================================
# 6. 导出文件
# ==========================================
plt.tight_layout()
plt.savefig('violin_plot_nature.pdf', bbox_inches='tight') # 导出供 AI 编辑的矢量图
plt.savefig('violin_plot_nature.png', dpi=300, bbox_inches='tight') # 导出高清预览图
plt.close()

print("Nature-style Violin Plot has been successfully saved!")