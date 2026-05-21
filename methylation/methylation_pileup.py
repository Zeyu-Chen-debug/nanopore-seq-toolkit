#!/usr/bin/env python3
"""Run modkit pileup on a BAM and export methylation calls to Excel.

Consolidates the per-modification scripts: choose the modified base with
``--mod-type`` (e.g. 5mC, 6mA), filter by a minimum percent-modified with
``--min-percent`` (0 keeps every site), and optionally reuse an existing
bgzipped BED with ``--skip-pileup``.

Example
-------
    python methylation_pileup.py \
        --bam SAMPLE.pass.bam \
        --reference b16_ref.fa \
        --mod-type 5mC \
        --out-dir b16c_meth/python_output \
        --min-percent 1
"""
import argparse
import subprocess
from pathlib import Path

import pandas as pd

BED_COLUMNS = [
    "chrom", "start_position", "end_position", "modified_base", "score",
    "strand", "start_position_2", "end_position_2", "color", "Nvalid_cov",
    "percent_modified", "Nmod", "Ncanonical", "Nother_mod", "Ndelete",
    "Nfail", "Ndiff", "Nnocall",
]


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-b", "--bam", required=True, type=Path,
                   help="Input modkit-compatible BAM with modification tags.")
    p.add_argument("-r", "--reference", required=True, type=Path,
                   help="Reference FASTA used by modkit pileup.")
    p.add_argument("-m", "--mod-type", default="5mC",
                   help="Modified base to call, e.g. 5mC or 6mA (default: 5mC).")
    p.add_argument("--out-dir", type=Path, default=None,
                   help="Output directory (created if needed). Used to derive default paths.")
    p.add_argument("--bed", type=Path, default=None,
                   help="Path for the bgzipped BED output (default: <out-dir>/SAMPLE.pass.bed.gz).")
    p.add_argument("--output-xlsx", type=Path, default=None,
                   help="Excel output path (default: <out-dir>/<mod-type>_methylation_sites.xlsx).")
    p.add_argument("--min-percent", type=float, default=1.0,
                   help="Keep sites with percent_modified greater than this (0 keeps all; default: 1).")
    p.add_argument("--modkit", default="modkit", help="Path to the modkit executable.")
    p.add_argument("--skip-pileup", action="store_true",
                   help="Reuse an existing BED (set with --bed) instead of running modkit.")
    return p


def resolve_paths(args):
    out_dir = args.out_dir or (args.bam.parent)
    out_dir.mkdir(parents=True, exist_ok=True)
    bed = args.bed or (out_dir / "SAMPLE.pass.bed.gz")
    xlsx = args.output_xlsx or (out_dir / f"{args.mod_type}_methylation_sites.xlsx")
    log = out_dir / "log.txt"
    return out_dir, bed, xlsx, log


def main(argv=None) -> None:
    args = build_parser().parse_args(argv)
    out_dir, bed, xlsx, log = resolve_paths(args)

    if not args.skip_pileup:
        cmd = [args.modkit, "pileup", str(args.bam), str(bed),
               "--modified-bases", args.mod_type,
               "--reference", str(args.reference),
               "--log", str(log), "--bgzf"]
        print("Running modkit pileup...")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.stdout:
            print("STDOUT:\n", result.stdout)
        if result.stderr:
            print("STDERR:\n", result.stderr)
        result.check_returncode()
        print(f"BED written to: {bed}")

    df = pd.read_csv(bed, sep="\t", names=BED_COLUMNS, compression="gzip")
    df["percent_modified"] = pd.to_numeric(df["percent_modified"], errors="coerce")
    if args.min_percent > 0:
        df = df[df["percent_modified"] > args.min_percent].copy()

    df.to_excel(xlsx, index=False)
    print(f"Wrote {len(df)} rows to: {xlsx}")


if __name__ == "__main__":
    main()
