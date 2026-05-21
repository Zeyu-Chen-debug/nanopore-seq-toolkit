#!/usr/bin/env python3
"""Print the mean error rate per group from an alignment summary spreadsheet.

The input Excel file is expected to contain an ``error_rate`` column and a
grouping column (default: ``source``). The mean error rate is printed for
each group value.

Example
-------
    python average_error_rate.py --input barcode0708.xlsx --group-col source
"""
import argparse
from pathlib import Path

import pandas as pd


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", required=True, type=Path,
                   help="Excel file with an 'error_rate' column.")
    p.add_argument("-g", "--group-col", default="source",
                   help="Column to group by (default: source).")
    p.add_argument("--value-col", default="error_rate",
                   help="Numeric column to average (default: error_rate).")
    return p


def main(argv=None) -> None:
    args = build_parser().parse_args(argv)
    df = pd.read_excel(args.input)
    means = df.groupby(args.group_col)[args.value_col].mean()
    for group, value in means.items():
        print(f"Average {args.value_col} for {group}: {value}")


if __name__ == "__main__":
    main()
