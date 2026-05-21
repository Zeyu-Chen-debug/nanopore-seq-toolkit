#!/usr/bin/env python3
"""Summarise read quality across one or more samples.

For each sample, reads are aligned to a reference with minimap2; a read is
counted as "good" if its alignment span covers at least a chosen fraction of
the target. The tool produces:

  * a 3-panel figure (good-read %, read-length distribution, error rate), and
  * an Excel table of the >=threshold-coverage alignments.

Samples are supplied as LABEL=FASTQ_DIR pairs, e.g.

    python good_reads_summary.py \
        --reference ref.fa \
        --sample FK=path/to/barcode05 \
        --sample FKS=path/to/barcode06 \
        --out-prefix results/N4_3h
"""
import argparse
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
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
    p.add_argument("-r", "--reference", required=True, type=Path,
                   help="Reference FASTA passed to minimap2.")
    p.add_argument("-s", "--sample", required=True, action="append",
                   metavar="LABEL=FASTQ_DIR",
                   help="A sample as LABEL=FASTQ_DIR. Repeat for more samples.")
    p.add_argument("-o", "--out-prefix", required=True, type=Path,
                   help="Output prefix; writes <prefix>.png and <prefix>.xlsx.")
    p.add_argument("--threshold", type=float, default=0.99,
                   help="Min fraction of target covered for a 'good' read (default: 0.99).")
    p.add_argument("--preset", default="map-ont",
                   help="minimap2 -x preset (default: map-ont).")
    p.add_argument("--glob", default="*.fastq.gz",
                   help="Glob used to select FASTQ files (default: *.fastq.gz).")
    p.add_argument("--minimap2", default="minimap2",
                   help="Path to the minimap2 executable (default: minimap2).")
    p.add_argument("--dpi", type=int, default=150, help="Figure DPI (default: 150).")
    return p


def parse_samples(pairs):
    samples = []
    for item in pairs:
        if "=" not in item:
            raise SystemExit(f"--sample must be LABEL=FASTQ_DIR, got: {item}")
        label, path = item.split("=", 1)
        samples.append((label, Path(path)))
    return samples


def run_minimap2(minimap2, preset, reference, fastq_dir, glob, paf_path):
    paf_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = (f'{minimap2} -x {preset} "{reference}" '
           f'{fastq_dir}/{glob} > "{paf_path}"')
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


def good_read_fraction(paf_path: Path, threshold: float):
    all_reads, good_reads = set(), set()
    with open(paf_path) as f:
        for line in f:
            fields = line.rstrip().split("\t")
            read_id = fields[0]
            tlen, tstart, tend = int(fields[6]), int(fields[7]), int(fields[8])
            all_reads.add(read_id)
            if abs(tend - tstart) >= threshold * tlen:
                good_reads.add(read_id)
    total = len(all_reads)
    pct = (len(good_reads) / total * 100) if total else 0.0
    return total, len(good_reads), pct


def main(argv=None) -> None:
    args = build_parser().parse_args(argv)
    samples = parse_samples(args.sample)
    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)

    labels, percentages, frames = [], [], []
    for label, fastq_dir in samples:
        paf_path = args.out_prefix.with_name(f"{args.out_prefix.name}_{label}.paf")
        run_minimap2(args.minimap2, args.preset, args.reference,
                     fastq_dir, args.glob, paf_path)
        total, good, pct = good_read_fraction(paf_path, args.threshold)
        print(f"{label}: {good}/{total} reads >= {args.threshold:.0%} coverage ({pct:.2f}%)")

        df = pd.read_csv(paf_path, sep="\t", header=None, usecols=range(12))
        df.columns = PAF_COLUMNS
        df["source"] = label
        df["alignment_span"] = df["target_end"] - df["target_start"]
        df["identity"] = df["residue_matches"] / df["alignment_block_length"]
        df["error_rate"] = 1 - df["identity"]
        df["cover_threshold"] = df["alignment_span"] / df["target_length"] >= args.threshold
        frames.append(df[df["cover_threshold"]])
        labels.append(label)
        percentages.append(pct)

    combined = pd.concat(frames, ignore_index=True)
    xlsx_path = args.out_prefix.with_suffix(".xlsx")
    combined.to_excel(xlsx_path, index=False)
    print(f"Excel table written to: {xlsx_path}")

    # ---------- figure ----------
    plt.style.use("seaborn-v0_8-whitegrid")
    palette = sns.color_palette("muted", n_colors=len(labels))
    color_of = dict(zip(labels, palette))
    fig, axes = plt.subplots(1, 3, figsize=(18, 5), dpi=args.dpi)

    # panel 1: good-read percentage
    ax = axes[0]
    bars = ax.bar(labels, percentages, color=[color_of[s] for s in labels])
    for bar, pct in zip(bars, percentages):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                f"{pct:.2f}%", ha="center", va="bottom", fontsize=10)
    ax.set_ylabel("Reads (%)")
    ax.set_title(f">= {args.threshold:.0%} insert coverage")
    ax.set_ylim(0, max(percentages) * 1.3 if percentages else 1)

    # panel 2: read-length density
    ax = axes[1]
    length_df = combined.drop_duplicates(subset=["source", "read_id"])
    for label in labels:
        subset = length_df[length_df["source"] == label]["read_length"]
        sns.kdeplot(subset, label=label, fill=True, alpha=0.3,
                    ax=ax, color=color_of[label])
    ax.set_xscale("log")
    ax.set_xlabel("Read length (bp)")
    ax.set_ylabel("Density")
    ax.set_title("Read length distribution")
    ax.legend(frameon=False)

    # panel 3: error-rate boxplot
    ax = axes[2]
    data = [combined[combined["source"] == s]["error_rate"].dropna() for s in labels]
    bp = ax.boxplot(data, labels=labels, showfliers=False, patch_artist=True,
                    widths=0.5, medianprops=dict(color="black", linewidth=2))
    for patch, label in zip(bp["boxes"], labels):
        patch.set_facecolor(color_of[label])
        patch.set_alpha(0.4)
    for i, label in enumerate(labels, start=1):
        y = combined[combined["source"] == label]["error_rate"].dropna().to_numpy()
        x = np.random.normal(i, 0.06, size=len(y))
        ax.scatter(x, y, s=10, alpha=0.3, color=color_of[label])
    ax.set_title("Error rate")
    ax.set_ylabel("Error rate")
    ax.set_ylim(0, 1)
    ax.grid(axis="y", linestyle="--", alpha=0.3)

    plt.tight_layout()
    png_path = args.out_prefix.with_suffix(".png")
    plt.savefig(png_path, dpi=args.dpi)
    print(f"Figure written to: {png_path}")


if __name__ == "__main__":
    main()
