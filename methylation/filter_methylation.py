#!/usr/bin/env python3
"""Filter a methylation Excel table by a minimum percent-modified value.

Example
-------
    python filter_methylation.py \
        --input barcode11_output.xlsx \
        --output high_methylation_sites.xlsx \
        --min-percent 10
"""
import argparse
from pathlib import Path

import pandas as pd


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", required=True, type=Path,
                   help="Excel file containing a percent-modified column.")
    p.add_argument("-o", "--output", required=True, type=Path,
                   help="Excel file to write the filtered rows to.")
    p.add_argument("--min-percent", type=float, default=10.0,
                   help="Keep rows with value strictly greater than this (default: 10).")
    p.add_argument("--column", default="percent modified",
                   help="Column to filter on (default: 'percent modified').")
    return p


def main(argv=None) -> None:
    args = build_parser().parse_args(argv)
    df = pd.read_excel(args.input)
    if args.column not in df.columns:
        raise SystemExit(f"Column '{args.column}' not found. Available: {list(df.columns)}")
    filtered = df[df[args.column] > args.min_percent]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    filtered.to_excel(args.output, index=False)
    print(f"Saved {len(filtered)} rows to {args.output}")


if __name__ == "__main__":
    main()
