#!/usr/bin/env python3
"""Plot per-base coverage depth across a reference for one sample.

Aligns reads with minimap2 to produce a PAF, accumulates depth over the target
positions, and writes a coverage line plot.

Example
-------
    python coverage_plot.py \
        --reference POC1610_ref.fa \
        --fastq-dir path/to/barcode11 \
        --output-paf barcode11.paf \
        --output barcode11_coverage.png \
        --title "barcode11 coverage across POC1610"
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import subprocess


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-r", "--reference", required=True, type=Path,
                   help="Reference FASTA passed to minimap2.")
    p.add_argument("-f", "--fastq-dir", required=True, type=Path,
                   help="Directory containing FASTQ reads.")
    p.add_argument("--output-paf", required=True, type=Path,
                   help="Path to write (or read with --skip-align) the PAF.")
    p.add_argument("-o", "--output", type=Path, default=None,
                   help="Output image path. If omitted, the figure is shown.")
    p.add_argument("--preset", default="map-ont",
                   help="minimap2 -x preset (default: map-ont).")
    p.add_argument("--glob", default="*.fastq.gz",
                   help="Glob used to select FASTQ files (default: *.fastq.gz).")
    p.add_argument("--minimap2", default="minimap2",
                   help="Path to the minimap2 executable (default: minimap2).")
    p.add_argument("--skip-align", action="store_true",
                   help="Use an existing PAF instead of re-running minimap2.")
    p.add_argument("--title", default="Coverage across reference", help="Plot title.")
    p.add_argument("--dpi", type=int, default=300, help="Figure DPI (default: 300).")
    return p


def main(argv=None) -> None:
    args = build_parser().parse_args(argv)

    if not args.skip_align:
        args.output_paf.parent.mkdir(parents=True, exist_ok=True)
        cmd = (f'{args.minimap2} -x {args.preset} "{args.reference}" '
               f'{args.fastq_dir}/{args.glob} > "{args.output_paf}"')
        subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")

    df = pd.read_csv(args.output_paf, sep="\t", header=None)
    ref_len = int(df[6].iloc[0])
    coverage = np.zeros(ref_len)
    for _, row in df.iterrows():
        coverage[int(row[7]):int(row[8])] += 1

    plt.figure(figsize=(10, 4))
    plt.plot(coverage)
    plt.xlabel("Reference position")
    plt.ylabel("Coverage")
    plt.title(args.title)
    plt.tight_layout()

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(args.output, dpi=args.dpi)
        print(f"Figure written to: {args.output}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
