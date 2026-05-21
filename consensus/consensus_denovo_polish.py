#!/usr/bin/env python3
"""Reference-free consensus by seeding on the longest read and polishing.

Picks the longest read as a seed, then iteratively aligns all reads back to it
(minimap2 + samtools) and applies bcftools variant calls for a configurable
number of polishing rounds.

Requires samtools and bcftools on PATH, plus a minimap2 executable.

Example
-------
    python consensus_denovo_polish.py \
        --fastq-dir path/to/barcode10 \
        --out-dir barcode10_consensus_denovo \
        --name barcode10 \
        --minimap2 ./minimap2/minimap2 \
        --rounds 2
"""
import argparse
import gzip
import shutil
import subprocess
from pathlib import Path
from typing import Iterator, Tuple


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-f", "--fastq-dir", required=True, type=Path,
                   help="Directory of FASTQ reads.")
    p.add_argument("-d", "--out-dir", required=True, type=Path,
                   help="Output directory (created if needed).")
    p.add_argument("-n", "--name", default="sample",
                   help="Basename for the final consensus (default: sample).")
    p.add_argument("--rounds", type=int, default=2,
                   help="Number of polishing rounds (default: 2).")
    p.add_argument("--preset", default="map-ont", help="minimap2 -ax preset (default: map-ont).")
    p.add_argument("--minimap2", default="minimap2",
                   help="Path to the minimap2 executable (default: minimap2).")
    return p


def open_maybe_gz(path: Path):
    return gzip.open(path, "rt") if path.name.endswith(".gz") else open(path, "rt")


def iterate_reads(path: Path) -> Iterator[Tuple[str, str]]:
    with open_maybe_gz(path) as fh:
        first = fh.readline()
        if not first:
            return
        fh.seek(0)
        if first.startswith(">"):
            name, parts = None, []
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if name is not None:
                        yield name, "".join(parts)
                    name, parts = line[1:].split()[0], []
                else:
                    parts.append(line)
            if name is not None:
                yield name, "".join(parts)
        elif first.startswith("@"):
            while True:
                name = fh.readline()
                if not name:
                    break
                seq = fh.readline()
                fh.readline()  # plus
                qual = fh.readline()
                if not qual:
                    break
                yield name[1:].strip().split()[0], seq.strip()
        else:
            raise ValueError(f"Unrecognised file format: {path}")


def write_fasta(path: Path, name: str, seq: str) -> None:
    with open(path, "w") as out:
        out.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            out.write(seq[i:i + 80] + "\n")


def collect_fastqs(fastq_dir: Path):
    files = sorted(
        list(fastq_dir.glob("*.fastq")) + list(fastq_dir.glob("*.fq")) +
        list(fastq_dir.glob("*.fastq.gz")) + list(fastq_dir.glob("*.fq.gz")))
    if not files:
        raise FileNotFoundError(f"No FASTQ files found in {fastq_dir}")
    return files


def pipe(cmd1, cmd2, err_label):
    p1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(cmd2, stdin=p1.stdout)
    p1.stdout.close()
    r2, r1 = p2.wait(), p1.wait()
    if r1 != 0 or r2 != 0:
        raise RuntimeError(f"{err_label} failed (codes {r1}, {r2})")


def main(argv=None) -> None:
    args = build_parser().parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    for tool in ("samtools", "bcftools"):
        if shutil.which(tool) is None:
            raise SystemExit(f"Required tool not found on PATH: {tool}")
    if shutil.which(args.minimap2) is None and not Path(args.minimap2).exists():
        raise SystemExit(f"minimap2 not found: {args.minimap2}")

    fastqs = collect_fastqs(args.fastq_dir)
    print(f"Found {len(fastqs)} FASTQ files")

    seed = args.out_dir / "seed_longest_read.fasta"
    bam = args.out_dir / "reads_to_seed.sorted.bam"
    vcf = args.out_dir / "variants.vcf.gz"
    final = args.out_dir / f"{args.name}_consensus.fasta"

    longest_name, longest_seq = None, ""
    for fq in fastqs:
        for name, seq in iterate_reads(fq):
            if len(seq) > len(longest_seq):
                longest_name, longest_seq = name, seq
    if not longest_seq:
        raise RuntimeError("No sequences were read from the FASTQ files.")
    print(f"Longest read: {longest_name} ({len(longest_seq)} bp)")
    write_fasta(seed, longest_name, longest_seq)

    current_ref = seed
    for rnd in range(1, args.rounds + 1):
        print(f"\n=== Polishing round {rnd} ===")
        consensus_path = args.out_dir / f"consensus_round{rnd}.fasta"
        pipe([args.minimap2, "-ax", args.preset, str(current_ref), *map(str, fastqs)],
             ["samtools", "sort", "-o", str(bam)], f"minimap2|sort round {rnd}")
        subprocess.run(["samtools", "index", str(bam)], check=True)
        pipe(["bcftools", "mpileup", "-Ou", "-f", str(current_ref), str(bam)],
             ["bcftools", "call", "-mv", "-Oz", "-o", str(vcf)], f"bcftools round {rnd}")
        subprocess.run(["bcftools", "index", str(vcf)], check=True)
        with open(consensus_path, "w") as out_fa:
            subprocess.run(["bcftools", "consensus", "-f", str(current_ref), str(vcf)],
                           check=True, stdout=out_fa)
        print(f"Consensus written to: {consensus_path}")
        current_ref = consensus_path

    shutil.copyfile(current_ref, final)
    print(f"\nFinal consensus FASTA written to: {final}")


if __name__ == "__main__":
    main()
