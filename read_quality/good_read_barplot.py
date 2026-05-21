#!/usr/bin/env python3
"""Grouped bar chart of a per-replicate metric, faceted by template.

Produces one subplot per template value, with bars grouped by treatment and
coloured by time point. Treatment, time and template values are read from the
data rather than hard-coded.

Example
-------
    python good_read_barplot.py --input cov_plt.xlsx --output good_reads.png
"""
import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", required=True, type=Path,
                   help="Excel file with facet, treatment, time and value columns.")
    p.add_argument("-o", "--output", type=Path, default=None,
                   help="Output image path. If omitted, the figure is shown.")
    p.add_argument("--facet-col", default="Plasmid",
                   help="Column to facet subplots by (default: Plasmid).")
    p.add_argument("--treatment-col", default="Trt",
                   help="Column for the x-axis groups (default: Trt).")
    p.add_argument("--time-col", default="Trt_time",
                   help="Column for the grouped bars (default: Trt_time).")
    p.add_argument("--value-col", default="pct_good_read",
                   help="Numeric column to plot (default: pct_good_read).")
    p.add_argument("--dpi", type=int, default=300, help="Figure DPI (default: 300).")
    return p


def main(argv=None) -> None:
    args = build_parser().parse_args(argv)
    df = pd.read_excel(args.input)

    facets = list(df[args.facet_col].dropna().unique())
    treatments = list(df[args.treatment_col].dropna().unique())
    times = list(df[args.time_col].dropna().unique())

    ncols = 2
    nrows = math.ceil(len(facets) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.5 * ncols, 4 * nrows),
                             sharey=True, squeeze=False)
    axes = axes.flatten()

    x = np.arange(len(treatments))
    width = 0.8 / max(len(times), 1)

    for ax, facet in zip(axes, facets):
        sub = df[df[args.facet_col] == facet]
        for j, time in enumerate(times):
            vals = []
            for trt in treatments:
                cell = sub[(sub[args.treatment_col] == trt) & (sub[args.time_col] == time)]
                vals.append(cell[args.value_col].iloc[0] if len(cell) else np.nan)
            offset = (j - (len(times) - 1) / 2) * width
            ax.bar(x + offset, vals, width=width, label=str(time))
        ax.set_title(str(facet))
        ax.set_xticks(x)
        ax.set_xticklabels([str(t) for t in treatments])
        ax.set_xlabel("Treatment")
        ax.grid(axis="y", alpha=0.3)

    for ax in axes[len(facets):]:
        ax.set_visible(False)
    axes[0].set_ylabel(args.value_col)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles[:len(times)], labels[:len(times)], title=args.time_col,
               loc="center left", bbox_to_anchor=(0.88, 0.5), frameon=False)
    fig.suptitle(f"{args.value_col} by {args.treatment_col}, {args.time_col} and {args.facet_col}",
                 y=0.98)
    fig.tight_layout(rect=[0, 0, 0.85, 0.95])

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(args.output, dpi=args.dpi)
        print(f"Figure written to: {args.output}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
