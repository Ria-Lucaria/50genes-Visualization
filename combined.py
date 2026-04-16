import os
from datetime import datetime
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

SUPPORTED_FORMATS = {"svg", "png", "pdf"}


def _load_and_prepare_data(csv_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV file not found: {csv_path}")
    if not csv_path.is_file():
        raise ValueError(f"CSV path is not a file: {csv_path}")

    df = pd.read_csv(csv_path)
    if df.empty:
        raise ValueError("CSV file is empty.")

    # 兼容历史列名格式
    if "Unnamed: 0" in df.columns:
        df = df.rename(columns={"Unnamed: 0": "Gene"})
    elif df.columns[0] == "":
        df.columns.values[0] = "Gene"
    required_cols = {"Gene", "NC", "dCas9"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"CSV must contain columns: {sorted(required_cols)}")

    df = df[["Gene", "NC", "dCas9"]].dropna().copy()
    if df.empty:
        raise ValueError("CSV has no valid rows after dropping NA values.")

    # 为了在线性轴上清晰显示大多数数据，去除掉 dCas9 组极端值
    upper_limit = df["dCas9"].quantile(0.95)
    df_trimmed = df[df["dCas9"] <= upper_limit].copy()
    if df_trimmed.empty:
        raise ValueError("No rows left after trimming top 5% dCas9 values.")

    df_melt = df_trimmed.melt(id_vars="Gene", value_vars=["NC", "dCas9"], var_name="Group", value_name="Reads")
    df_melt["Group"] = df_melt["Group"].replace({"dCas9": "WORF"})

    df_trimmed["Diff"] = df_trimmed["dCas9"] - df_trimmed["NC"]
    df_sorted = df_trimmed.sort_values(by="Diff", ascending=False).reset_index(drop=True)
    heat_df = df_sorted[["Gene", "NC", "dCas9"]].set_index("Gene")
    heat_df.columns = ["NC", "WORF"]

    return df_melt, heat_df


def render_combined_figure(
    csv_path: str,
    output_dir: str = "output",
    filename_prefix: str = "Fig6",
    output_format: str = "svg",
) -> str:
    output_format = output_format.lower()
    if output_format not in SUPPORTED_FORMATS:
        raise ValueError(f"Unsupported output format: {output_format}. Supported: {sorted(SUPPORTED_FORMATS)}")
    if not filename_prefix or not filename_prefix.strip():
        raise ValueError("filename_prefix must not be empty.")

    csv_path_obj = Path(csv_path).expanduser().resolve()
    output_dir_obj = Path(output_dir).expanduser().resolve()
    os.makedirs(output_dir_obj, exist_ok=True)

    if not os.access(output_dir_obj, os.W_OK):
        raise PermissionError(f"Output directory is not writable: {output_dir_obj}")

    df_melt, heat_df = _load_and_prepare_data(csv_path_obj)

    mpl.rcParams["font.family"] = "Arial"
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["axes.linewidth"] = 1.0

    fig, (ax1, ax2) = plt.subplots(
        1,
        2,
        figsize=(8, 5),
        dpi=300,
        gridspec_kw={"width_ratios": [2, 1], "wspace": 0.35},
    )

    colors = {"NC": "#A0CBE8", "WORF": "#FF9D9A"}
    sns.violinplot(
        x="Group",
        y="Reads",
        data=df_melt,
        palette=colors,
        inner=None,
        linewidth=0,
        alpha=0.6,
        ax=ax1,
        zorder=1,
    )
    sns.boxplot(
        x="Group",
        y="Reads",
        data=df_melt,
        width=0.15,
        boxprops={"facecolor": "none", "edgecolor": "black", "zorder": 2, "linewidth": 0.5},
        medianprops={"color": "black", "linewidth": 1, "zorder": 2},
        whiskerprops={"linewidth": 0.5, "color": "black"},
        capprops={"linewidth": 0.5, "color": "black"},
        showfliers=False,
        ax=ax1,
    )
    sns.stripplot(
        x="Group",
        y="Reads",
        data=df_melt,
        palette=colors,
        size=2.5,
        jitter=0.15,
        alpha=0.8,
        ax=ax1,
        zorder=3,
        edgecolor="none",
    )

    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.set_xlabel("")
    ax1.set_ylabel("Reads", fontsize=12)
    ax1.tick_params(axis="both", which="major", labelsize=11)

    cmap = sns.color_palette("rocket_r", as_cmap=True)
    sns.heatmap(
        heat_df,
        cmap=cmap,
        ax=ax2,
        yticklabels=False,
        cbar_kws={"label": "Reads", "shrink": 0.7, "aspect": 30},
        linewidths=0,
    )

    ax2.set_ylabel("Target Regions (Ranked by Absolute Enrichment)", fontsize=12)
    ax2.tick_params(axis="x", which="major", labelsize=11, rotation=0)
    cbar = ax2.collections[0].colorbar
    cbar.ax.tick_params(labelsize=9)
    cbar.set_label("Reads", size=11)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_path = output_dir_obj / f"{filename_prefix}+{timestamp}.{output_format}"
    fig.savefig(output_path, bbox_inches="tight", format=output_format)
    plt.close(fig)

    return str(output_path)


def main() -> None:
    default_csv = str(Path.cwd() / "all_reads_compare.csv")
    output_file = render_combined_figure(default_csv)
    print(f"已剔除Top 5%的极端异常值，图表已保存为 {output_file}")


if __name__ == "__main__":
    main()
