#!/usr/bin/env python3
"""Draw each read's alignment span along the reference.

Useful for visually inspecting where reads start and end on the target (for
example to spot linearisation/cut sites). Reads an existing PAF, or generates
one first with minimap2.

Example
-------
    # from an existing PAF
    python alignment_spans.py --paf barcode17.paf --output spans.png

    # generate the PAF first
    python alignment_spans.py --reference b17_ref.fa \
        --fastq-dir path/to/barcode17 --output-paf barcode17.paf
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import subprocess


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--paf", type=Path, default=None,
                   help="Existing PAF to plot. If omitted, --reference and "
                        "--fastq-dir are used to generate one.")
    p.add_argument("-r", "--reference", type=Path, default=None,
                   help="Reference FASTA (needed when generating a PAF).")
    p.add_argument("-f", "--fastq-dir", type=Path, default=None,
                   help="FASTQ directory (needed when generating a PAF).")
    p.add_argument("--output-paf", type=Path, default=None,
                   help="Where to write a generated PAF.")
    p.add_argument("-o", "--output", type=Path, default=None,
                   help="Output image path. If omitted, the figure is shown.")
    p.add_argument("--preset", default="map-ont", help="minimap2 -x preset (default: map-ont).")
    p.add_argument("--glob", default="*.fastq.gz", help="FASTQ glob (default: *.fastq.gz).")
    p.add_argument("--minimap2", default="minimap2", help="minimap2 executable.")
    p.add_argument("--title", default="Read alignments along reference", help="Plot title.")
    p.add_argument("--dpi", type=int, default=300, help="Figure DPI (default: 300).")
    return p


def main(argv=None) -> None:
    args = build_parser().parse_args(argv)

    paf_path = args.paf
    if paf_path is None:
        if not (args.reference and args.fastq_dir and args.output_paf):
            raise SystemExit("Provide --paf, or all of --reference, --fastq-dir "
                             "and --output-paf to generate one.")
        args.output_paf.parent.mkdir(parents=True, exist_ok=True)
        cmd = (f'{args.minimap2} -x {args.preset} "{args.reference}" '
               f'{args.fastq_dir}/{args.glob} > "{args.output_paf}"')
        subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
        paf_path = args.output_paf

    df = pd.read_csv(paf_path, sep="\t", header=None)
    target_start = pd.to_numeric(df[7]).to_numpy()
    target_end = pd.to_numeric(df[8]).to_numpy()

    plt.figure(figsize=(10, 6))
    for i in range(len(df)):
        plt.plot([target_start[i], target_end[i]], [i, i], linewidth=0.5)
    plt.xlabel("Reference (target) position")
    plt.ylabel("Read index")
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
