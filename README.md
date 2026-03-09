# 50genes-Visualization

50个基因 reads 差异可视化。
注意：all_reads_comapre.csv是实验室数据，不在这里提供。本repo只分享绘图方法。

## 更新日志

### 2026-03-09

1. 图版布局调整
- `combined.py` 中 Panel 宽度比例已调整为 `A:B = 2:1`。
- 说明：当前代码使用 `gridspec_kw={'width_ratios': [2, 1], 'wspace': 0.35}`。

2. 分组名称显示调整
- 图中原 `dCas9` 显示名改为 `WORF`（仅改图中显示，不改原始数据列名）。
- 小提琴图与热图列名显示已同步为 `NC` 和 `WORF`。

3. 导出格式与命名规则调整
- 导出格式从 `PNG/PDF` 改为仅导出 `SVG`。
- 输出目录改为项目下 `output/`（不存在时自动创建）。
- 文件名改为 `Fig6+时间戳.svg`，例如 `Fig6+20260309_133353.svg`。

4. Panel 标签处理
- 已移除自动添加的 `A`、`B` 标签，方便后续手动排版。

5. 脚本运行与环境验证
- 已在 `base` conda 环境完成运行测试。
- 为满足运行依赖，已安装 `seaborn`。
- 脚本可正常执行并在 `output/` 目录生成目标 SVG 文件。

6. 当前已知提示（不影响出图）
- `seaborn` 关于 `palette` 用法存在 FutureWarning。
- 系统未安装 `Arial` 时会出现字体回退提示（`findfont: Font family 'Arial' not found.`）。
