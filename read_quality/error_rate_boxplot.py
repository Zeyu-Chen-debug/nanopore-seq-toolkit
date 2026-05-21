#!/usr/bin/env python3
"""Box-and-strip plot of per-read error rate across alignment groups.

Each group is an existing PAF file supplied as LABEL=PAF. Reads are filtered
to those covering at least ``--threshold`` of the target before plotting.

Example
-------
    python error_rate_boxplot.py \
        --paf "FK (30 min)=N1_FK_30min.paf" \
        --paf "FK (3 h)=N1_FK_3h.paf" \
        --paf "FKS (30 min)=N1_FKS_30min.paf" \
        --paf "FKS (3 h)=N1_FKS_3h.paf" \
        --title "Error rate N1" --output N1_error.png
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
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
                   help="Output image path (e.g. error.png). Omit to show interactively.")
    p.add_argument("--threshold", type=float, default=0.99,
                   help="Min fraction of target covered to keep a read (default: 0.99).")
    p.add_argument("--title", default="Error rate", help="Plot title.")
    p.add_argument("--dpi", type=int, default=300, help="Figure DPI (default: 300).")
    p.add_argument("--show", action="store_true",
                   help="Display the figure instead of (or as well as) saving.")
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
        df["identity"] = df["residue_matches"] / df["alignment_block_length"]
        df["error_rate"] = 1 - df["identity"]
        df["cover_ok"] = df["alignment_span"] / df["target_length"] >= threshold
        frames.append(df[df["cover_ok"]])
    return pd.concat(frames, ignore_index=True)


def main(argv=None) -> None:
    args = build_parser().parse_args(argv)
    df = load_groups(args.paf, args.threshold)
    order = [item.split("=", 1)[0] for item in args.paf]

    plt.figure(figsize=(8, 6))
    sns.boxplot(data=df, x="group_label", y="error_rate", order=order,
                width=0.6, showfliers=False)
    sns.stripplot(data=df, x="group_label", y="error_rate", order=order,
                  color="black", alpha=0.4, jitter=0.25, size=3)
    plt.title(args.title, fontsize=14)
    plt.xlabel("")
    plt.ylabel("Error rate")
    plt.ylim(0, 1)
    plt.xticks(rotation=20)
    sns.despine()
    plt.tight_layout()

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(args.output, dpi=args.dpi)
        print(f"Figure written to: {args.output}")
    if args.show or not args.output:
        plt.show()


if __name__ == "__main__":
    main()
