#!/usr/bin/env python3
"""Scatter of good-read percentage against template length, with a trendline.

Reads a spreadsheet with one row per replicate and plots ``pct_good_read``
versus a length column, colouring points by a grouping column and overlaying a
linear fit across all points.

Example
-------
    python good_read_vs_length.py --input cov_plt.xlsx --output good_vs_len.png
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", required=True, type=Path,
                   help="Excel file with length, group and value columns.")
    p.add_argument("-o", "--output", type=Path, default=None,
                   help="Output image path. If omitted, the figure is shown.")
    p.add_argument("--length-col", default="Plasmid_length",
                   help="Column with the template length (default: Plasmid_length).")
    p.add_argument("--value-col", default="pct_good_read",
                   help="Column with the value to plot (default: pct_good_read).")
    p.add_argument("--group-col", default="Plasmid",
                   help="Column used to colour points (default: Plasmid).")
    p.add_argument("--jitter", type=float, default=10.0,
                   help="Std-dev of horizontal jitter added to points (default: 10).")
    p.add_argument("--dpi", type=int, default=300, help="Figure DPI (default: 300).")
    return p


def main(argv=None) -> None:
    args = build_parser().parse_args(argv)
    df = pd.read_excel(args.input)

    plt.figure(figsize=(6, 4))
    for group, sub in df.groupby(args.group_col):
        x = sub[args.length_col].iloc[0]
        jitter = np.random.normal(0, args.jitter, size=len(sub))
        plt.scatter(x + jitter, sub[args.value_col], label=group)

    x_all = df[args.length_col].values
    y_all = df[args.value_col].values
    m, b = np.polyfit(x_all, y_all, 1)
    x_line = np.linspace(min(x_all), max(x_all), 100)
    plt.plot(x_line, m * x_line + b, linestyle="--", linewidth=2)

    plt.xlabel("Template length")
    plt.ylabel("Percent of good reads")
    plt.title(f"{args.value_col} vs {args.length_col}")
    plt.legend(title=args.group_col)
    plt.tight_layout()

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(args.output, dpi=args.dpi)
        print(f"Figure written to: {args.output}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
