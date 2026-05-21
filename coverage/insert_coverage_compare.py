#!/usr/bin/env python3
"""Compare insert coverage across several samples, with optional subsampling.

For each sample a random subset of reads is drawn (for a fair comparison at
equal depth), aligned to the reference, and the per-base coverage is overlaid
on a single plot.

Example
-------
    python insert_coverage_compare.py \
        --reference N4T_w_adaptor.fa \
        --sample "FK 30min=path/to/barcode14" \
        --sample "FK 3h=path/to/barcode07" \
        --sample "FKS 30min=path/to/barcode15" \
        --sample "FKS 3h=path/to/barcode08" \
        --n-reads 2000 --output Coverage_N4.png
"""
import argparse
import gzip
import random
import subprocess
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-r", "--reference", required=True, type=Path,
                   help="Reference FASTA passed to minimap2.")
    p.add_argument("-s", "--sample", required=True, action="append",
                   metavar="LABEL=FASTQ_DIR",
                   help="A sample as LABEL=FASTQ_DIR. Repeat for more samples.")
    p.add_argument("--paf-dir", type=Path, default=Path("."),
                   help="Directory to write intermediate PAF files (default: .).")
    p.add_argument("-o", "--output", type=Path, default=None,
                   help="Output image path. If omitted, the figure is shown.")
    p.add_argument("--n-reads", type=int, default=2000,
                   help="Reads to subsample per sample (default: 2000).")
    p.add_argument("--seed", type=int, default=42, help="Random seed (default: 42).")
    p.add_argument("--preset", default="map-ont",
                   help="minimap2 -x preset (default: map-ont).")
    p.add_argument("--glob", default="*.fastq.gz",
                   help="Glob used to select FASTQ files (default: *.fastq.gz).")
    p.add_argument("--minimap2", default="minimap2",
                   help="Path to the minimap2 executable (default: minimap2).")
    p.add_argument("--ymax", type=float, default=None, help="Optional y-axis upper limit.")
    p.add_argument("--title", default="Insert coverage", help="Plot title.")
    p.add_argument("--dpi", type=int, default=150, help="Figure DPI (default: 150).")
    return p


def fasta_length(path: Path) -> int:
    length = 0
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                length += len(line.strip())
    return length


def iter_fastq_records(path: Path):
    with gzip.open(path, "rt") as f:
        while True:
            h = f.readline()
            if not h:
                break
            yield h + f.readline() + f.readline() + f.readline()


def subsample(fastq_dir: Path, glob: str, n_reads: int, out_fastq: Path, rng: random.Random) -> int:
    records = []
    fastqs = sorted(fastq_dir.glob(glob))
    if not fastqs:
        raise FileNotFoundError(f"No FASTQ files in {fastq_dir}")
    for fq in fastqs:
        records.extend(iter_fastq_records(fq))
    selected = rng.sample(records, min(n_reads, len(records)))
    with open(out_fastq, "w") as out:
        out.writelines(selected)
    return len(selected)


def paf_to_coverage(paf_path: Path, target_len: int) -> np.ndarray:
    diff = np.zeros(target_len + 1, dtype=int)
    with open(paf_path) as f:
        for line in f:
            fields = line.split("\t")
            tstart, tend = int(fields[7]), int(fields[8])
            if tend > tstart:
                diff[tstart] += 1
                diff[tend] -= 1
    return diff.cumsum()[:-1]


def main(argv=None) -> None:
    args = build_parser().parse_args(argv)
    rng = random.Random(args.seed)
    args.paf_dir.mkdir(parents=True, exist_ok=True)
    target_len = fasta_length(args.reference)

    plt.figure(figsize=(10, 4))
    for item in args.sample:
        if "=" not in item:
            raise SystemExit(f"--sample must be LABEL=FASTQ_DIR, got: {item}")
        label, path = item.split("=", 1)
        fastq_dir = Path(path)
        paf_path = args.paf_dir / f"{label.replace(' ', '_')}.paf"

        with tempfile.TemporaryDirectory() as tmp:
            subset = Path(tmp) / "subset.fastq"
            n = subsample(fastq_dir, args.glob, args.n_reads, subset, rng)
            print(f"{label}: {n} reads used")
            cmd = [args.minimap2, "-x", args.preset, str(args.reference), str(subset)]
            with open(paf_path, "w") as out:
                result = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
            if result.returncode != 0:
                raise RuntimeError(result.stderr)

        coverage = paf_to_coverage(paf_path, target_len)
        plt.plot(np.arange(1, target_len + 1), coverage, label=label, linewidth=1.5)

    if args.ymax is not None:
        plt.ylim(0, args.ymax)
    plt.xlabel("Position on reference")
    plt.ylabel("Coverage depth")
    plt.title(args.title)
    plt.legend()
    plt.tight_layout()

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(args.output, dpi=args.dpi)
        print(f"Figure written to: {args.output}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
