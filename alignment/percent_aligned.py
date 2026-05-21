#!/usr/bin/env python3
"""Report the percentage of reads that aligned to a reference.

Runs minimap2 to produce a PAF, counts the total reads in the FASTQ
directory, counts the unique aligned reads in the PAF, and prints the
alignment rate.

Example
-------
    python percent_aligned.py \
        --reference ref.fa \
        --fastq-dir path/to/barcode08 \
        --output-paf barcode08.paf
"""
import argparse
import gzip
import subprocess
from pathlib import Path

import pandas as pd


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-r", "--reference", required=True, type=Path,
                   help="Reference FASTA passed to minimap2.")
    p.add_argument("-f", "--fastq-dir", required=True, type=Path,
                   help="Directory containing gzipped FASTQ reads.")
    p.add_argument("-o", "--output-paf", required=True, type=Path,
                   help="Path to write the PAF output.")
    p.add_argument("--preset", default="map-ont",
                   help="minimap2 -x preset (default: map-ont).")
    p.add_argument("--glob", default="*.fastq.gz",
                   help="Glob used to select FASTQ files (default: *.fastq.gz).")
    p.add_argument("--minimap2", default="minimap2",
                   help="Path to the minimap2 executable (default: minimap2).")
    p.add_argument("--skip-align", action="store_true",
                   help="Use an existing PAF instead of re-running minimap2.")
    return p


def count_reads(fastq_dir: Path, glob: str) -> int:
    total = 0
    for fq in fastq_dir.glob(glob):
        with gzip.open(fq, "rt") as f:
            total += sum(1 for _ in f) // 4
    return total


def main(argv=None) -> None:
    args = build_parser().parse_args(argv)

    if not args.skip_align:
        args.output_paf.parent.mkdir(parents=True, exist_ok=True)
        cmd = (f'{args.minimap2} -x {args.preset} "{args.reference}" '
               f'{args.fastq_dir}/{args.glob} > "{args.output_paf}"')
        subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")

    total_reads = count_reads(args.fastq_dir, args.glob)
    if total_reads == 0:
        raise SystemExit("ERROR: No reads found in FASTQ files.")

    paf_df = pd.read_csv(args.output_paf, sep="\t", header=None)
    aligned_reads = paf_df[0].nunique()
    pct = aligned_reads / total_reads * 100

    print("\n===== Alignment Summary =====")
    print(f"Total reads:        {total_reads:,}")
    print(f"Aligned reads:      {aligned_reads:,}")
    print(f"Percentage aligned: {pct:.2f}%")


if __name__ == "__main__":
    main()
