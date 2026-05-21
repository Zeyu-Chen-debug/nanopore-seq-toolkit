#!/usr/bin/env python3
"""Grid of box-and-strip plots summarising QC metrics by treatment and time.

Column names are normalised on load (spaces -> underscores, parentheses
removed, '%' -> 'pct'). Each requested metric is drawn as a boxplot split by a
treatment column and a hue (time) column.

Example
-------
    python summary_boxplots.py --input cov_plt.xlsx \
        --metrics N50kb mean_coverage pct_good_read av_error_rate \
        --output summary.png
"""
import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", required=True, type=Path,
                   help="Excel file of QC metrics, one row per replicate.")
    p.add_argument("-o", "--output", type=Path, default=None,
                   help="Output image path. If omitted, the figure is shown.")
    p.add_argument("-m", "--metrics", nargs="+",
                   default=["N50kb", "mean_coverage", "pct_good_read", "av_error_rate"],
                   help="Normalised metric column names to plot.")
    p.add_argument("--treatment-col", default="Trt",
                   help="Column for the x-axis (default: Trt).")
    p.add_argument("--hue-col", default="Trt_time",
                   help="Column for the hue split (default: Trt_time).")
    p.add_argument("--dpi", type=int, default=300, help="Figure DPI (default: 300).")
    return p


def normalise_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = (df.columns.str.strip()
                  .str.replace(" ", "_", regex=False)
                  .str.replace("(", "", regex=False)
                  .str.replace(")", "", regex=False)
                  .str.replace("%", "pct", regex=False))
    return df


def main(argv=None) -> None:
    args = build_parser().parse_args(argv)
    df = normalise_columns(pd.read_excel(args.input))
    df[args.treatment_col] = df[args.treatment_col].astype("category")
    df[args.hue_col] = df[args.hue_col].astype("category")

    missing = [m for m in args.metrics if m not in df.columns]
    if missing:
        raise SystemExit(f"Columns not found after normalisation: {missing}\n"
                         f"Available: {list(df.columns)}")

    ncols = 2
    nrows = math.ceil(len(args.metrics) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 5 * nrows), squeeze=False)
    axes = axes.flatten()

    for ax, metric in zip(axes, args.metrics):
        sns.boxplot(data=df, x=args.treatment_col, y=metric, hue=args.hue_col, ax=ax)
        sns.stripplot(data=df, x=args.treatment_col, y=metric, hue=args.hue_col,
                      dodge=True, color="black", alpha=0.6, ax=ax, legend=False)
        ax.set_title(metric)
        ax.set_xlabel("Treatment")
        ax.set_ylabel(metric)
        handles, labels = ax.get_legend_handles_labels()
        n = df[args.hue_col].nunique()
        ax.legend(handles[:n], labels[:n], title=args.hue_col)

    for ax in axes[len(args.metrics):]:
        ax.set_visible(False)

    plt.tight_layout()
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(args.output, dpi=args.dpi)
        print(f"Figure written to: {args.output}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
