#!/usr/bin/env python3
"""Align Nanopore reads to a reference with minimap2 and write a PAF file.

Optionally also writes an Excel summary of the alignments with derived
per-alignment metrics (alignment span, identity, mismatches, error rate).

Example
-------
    python generate_paf.py \
        --reference ref.fa \
        --fastq-dir path/to/barcode01 \
        --output-paf barcode01.paf \
        --to-xlsx
"""
import argparse
import subprocess
from pathlib import Path

import pandas as pd

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
    p.add_argument("-f", "--fastq-dir", required=True, type=Path,
                   help="Directory containing FASTQ reads.")
    p.add_argument("-o", "--output-paf", required=True, type=Path,
                   help="Path to write the PAF output.")
    p.add_argument("--preset", default="map-ont",
                   help="minimap2 -x preset (default: map-ont).")
    p.add_argument("--glob", default="*.fastq.gz",
                   help="Glob used to select FASTQ files (default: *.fastq.gz).")
    p.add_argument("--minimap2", default="minimap2",
                   help="Path to the minimap2 executable (default: minimap2).")
    p.add_argument("--to-xlsx", action="store_true",
                   help="Also write an .xlsx summary next to the PAF.")
    return p


def run_minimap2(minimap2: str, preset: str, reference: Path,
                 fastq_dir: Path, glob: str, output_paf: Path) -> None:
    output_paf.parent.mkdir(parents=True, exist_ok=True)
    cmd = (f'{minimap2} -x {preset} "{reference}" '
           f'{fastq_dir}/{glob} > "{output_paf}"')
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


def paf_to_dataframe(output_paf: Path) -> pd.DataFrame:
    df = pd.read_csv(output_paf, sep="\t", header=None, usecols=range(12))
    df.columns = PAF_COLUMNS
    df["alignment_span"] = df["target_end"] - df["target_start"]
    df["identity"] = df["residue_matches"] / df["alignment_block_length"]
    df["mismatches"] = df["alignment_block_length"] - df["residue_matches"]
    df["error_rate"] = 1 - df["identity"]
    return df


def main(argv=None) -> None:
    args = build_parser().parse_args(argv)
    run_minimap2(args.minimap2, args.preset, args.reference,
                 args.fastq_dir, args.glob, args.output_paf)
    print(f"PAF written to: {args.output_paf}")
    if args.to_xlsx:
        df = paf_to_dataframe(args.output_paf)
        xlsx_path = args.output_paf.with_suffix(".xlsx")
        df.to_excel(xlsx_path, index=False)
        print(f"Excel summary written to: {xlsx_path}")


if __name__ == "__main__":
    main()
