#!/usr/bin/env python3
"""Plot the distribution of mismatch positions along a reference.

Streams FASTQ reads through minimap2 with the ``--cs`` tag, keeps alignments
passing identity and coverage thresholds, parses the cs tag to locate
substitutions, and plots a per-position mismatch count.

Example
-------
    python mismatch_positions.py \
        --reference N1T_w_adaptor.fa \
        --reads-dir path/to/barcode01 \
        --output-paf output.paf \
        --output mismatches.png
"""
import argparse
import glob
import gzip
import os
import re
import subprocess
from collections import defaultdict

import matplotlib.pyplot as plt


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-r", "--reference", required=True,
                   help="Reference FASTA passed to minimap2.")
    p.add_argument("-f", "--reads-dir", required=True,
                   help="Directory of gzipped FASTQ reads.")
    p.add_argument("--output-paf", default="output.paf",
                   help="PAF path to write/read (default: output.paf).")
    p.add_argument("-o", "--output", default=None,
                   help="Output image path. If omitted, the figure is shown.")
    p.add_argument("--min-identity", type=float, default=0.99,
                   help="Minimum alignment identity to keep (default: 0.99).")
    p.add_argument("--min-coverage", type=float, default=0.99,
                   help="Minimum target coverage to keep (default: 0.99).")
    p.add_argument("--threads", type=int, default=8, help="minimap2 threads (default: 8).")
    p.add_argument("--preset", default="map-ont", help="minimap2 -x preset (default: map-ont).")
    p.add_argument("--glob", default="*.fastq.gz", help="FASTQ glob (default: *.fastq.gz).")
    p.add_argument("--minimap2", default="minimap2", help="minimap2 executable.")
    p.add_argument("--reuse-paf", action="store_true",
                   help="Reuse an existing PAF instead of regenerating it.")
    return p


def generate_paf(args) -> None:
    if os.path.exists(args.output_paf) and args.reuse_paf:
        print(f"[INFO] Reusing existing PAF: {args.output_paf}")
        return
    if os.path.exists(args.output_paf):
        os.remove(args.output_paf)

    fastq_files = glob.glob(os.path.join(args.reads_dir, args.glob))
    if not fastq_files:
        raise ValueError(f"No {args.glob} files found in {args.reads_dir}")
    print(f"[INFO] Found {len(fastq_files)} FASTQ files")

    cmd = [args.minimap2, "-t", str(args.threads), "--cs", "-x", args.preset,
           args.reference, "-"]
    with open(args.output_paf, "w") as out:
        proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=out, text=True)
        for fq in fastq_files:
            print(f"[INFO] Streaming {fq}")
            with gzip.open(fq, "rt") as f:
                for line in f:
                    proc.stdin.write(line)
        proc.stdin.close()
        proc.wait()
    print(f"[INFO] PAF generated: {args.output_paf}")


def parse_paf(paf_file: str, min_identity: float, min_coverage: float):
    mismatch_counts = defaultdict(int)
    total = passed = 0
    with open(paf_file) as f:
        for line in f:
            total += 1
            fields = line.strip().split("\t")
            tlen = int(fields[6])
            tstart, tend = int(fields[7]), int(fields[8])
            matches, aln_len = int(fields[9]), int(fields[10])
            identity = matches / aln_len if aln_len else 0
            coverage = (tend - tstart) / tlen if tlen else 0
            if identity < min_identity or coverage < min_coverage:
                continue
            passed += 1
            cs_tag = next((fld[5:] for fld in fields[12:] if fld.startswith("cs:Z:")), None)
            if cs_tag is None:
                continue
            t_pos = tstart
            for token in re.findall(r'(:\d+|\*[acgtn][acgtn]|[+-][acgtn]+)', cs_tag):
                if token.startswith(':'):
                    t_pos += int(token[1:])
                elif token.startswith('*'):
                    mismatch_counts[t_pos] += 1
                    t_pos += 1
                elif token.startswith('-'):
                    t_pos += len(token) - 1
                # '+' (insertion) does not advance the target position
    print(f"[INFO] Total alignments: {total}")
    print(f"[INFO] Passed filters:   {passed}")
    print(f"[INFO] Mismatch positions: {len(mismatch_counts)}")
    return mismatch_counts


def main(argv=None) -> None:
    args = build_parser().parse_args(argv)
    generate_paf(args)
    mismatch_counts = parse_paf(args.output_paf, args.min_identity, args.min_coverage)
    if not mismatch_counts:
        print("[WARNING] No mismatches found - nothing to plot.")
        return

    x = sorted(mismatch_counts)
    y = [mismatch_counts[pos] for pos in x]
    plt.figure(figsize=(12, 5))
    plt.scatter(x, y, s=10)
    plt.xlabel("Target position")
    plt.ylabel("Mismatch count")
    plt.title(f"Mismatch distribution (>= {args.min_identity:.0%} identity & coverage)")
    plt.grid(True)
    plt.tight_layout()

    if args.output:
        plt.savefig(args.output, dpi=300)
        print(f"Figure written to: {args.output}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
