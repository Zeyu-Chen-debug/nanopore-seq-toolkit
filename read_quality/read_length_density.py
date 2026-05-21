#!/usr/bin/env python3
"""Kernel-density plot of read-length distributions across alignment groups.

Each group is an existing PAF file supplied as LABEL=PAF. Reads are filtered
to those covering at least ``--threshold`` of the target before plotting, and
read length is shown on a log axis.

Example
-------
    python read_length_density.py \
        --paf "FK (30 min)=N4_FK_30min.paf" \
        --paf "FKS (3 h)=N4_FKS_3h.paf" \
        --title "Read length distribution N4" --output N4_rl.png
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import pandas as pd
import seaborn as sns

PAF_COLUMNS = [
    "read_id", "read_length", "read_start", "read_end", "strand",
    "target_name", "target_length", "target_start", "target_end",
    "residue_matches", "alignment_block_length", "mapping_quality",
]


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-p", "--paf", required=True, action="append",
                   metavar="LABEL=PAF",
                   help="A group as LABEL=PAF_PATH. Repeat for more groups.")
    p.add_argument("-o", "--output", required=True, type=Path,
                   help="Output image path (e.g. rl.png).")
    p.add_argument("--threshold", type=float, default=0.99,
                   help="Min fraction of target covered to keep a read (default: 0.99).")
    p.add_argument("--title", default="Read length distribution", help="Plot title.")
    p.add_argument("--bw-adjust", type=float, default=0.9,
                   help="KDE bandwidth adjustment (default: 0.9).")
    p.add_argument("--dpi", type=int, default=300, help="Figure DPI (default: 300).")
    p.add_argument("--show", action="store_true", help="Also display the figure.")
    return p


def load_groups(pairs, threshold):
    frames = []
    for item in pairs:
        if "=" not in item:
            raise SystemExit(f"--paf must be LABEL=PAF, got: {item}")
        label, path = item.split("=", 1)
        df = pd.read_csv(path, sep="\t", header=None, usecols=range(12))
        df.columns = PAF_COLUMNS
        df["group_label"] = label
        df["alignment_span"] = df["target_end"] - df["target_start"]
        df["cover_ok"] = df["alignment_span"] / df["target_length"] >= threshold
        frames.append(df[df["cover_ok"]])
    return pd.concat(frames, ignore_index=True)


def main(argv=None) -> None:
    args = build_parser().parse_args(argv)
    df = load_groups(args.paf, args.threshold)
    order = [item.split("=", 1)[0] for item in args.paf]

    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(9, 6))
    sns.kdeplot(data=df, x="read_length", hue="group_label", hue_order=order,
                fill=True, alpha=0.35, linewidth=1.5, bw_adjust=args.bw_adjust,
                common_norm=False, ax=ax)
    ax.set_xscale("log")
    ax.set_xlabel("Read length (bp)", fontsize=12)
    ax.set_ylabel("Density", fontsize=12)
    ax.set_title(args.title, fontsize=14)
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    legend = ax.get_legend()
    if legend is not None:
        handles = legend.legend_handles
        labels = [t.get_text() for t in legend.texts]
        legend.remove()
        fig.legend(handles, labels, title="Condition", loc="center right", frameon=True)
        plt.subplots_adjust(right=0.8)

    sns.despine()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(args.output, dpi=args.dpi)
    print(f"Figure written to: {args.output}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
