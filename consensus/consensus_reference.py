#!/usr/bin/env python3
"""Build a reference-guided consensus from Nanopore reads.

Aligns reads to a reference (minimap2), sorts/indexes (samtools), calls
variants (bcftools mpileup/call) and applies them to produce a consensus FASTA.

Requires minimap2, samtools and bcftools on PATH.

Example
-------
    python consensus_reference.py \
        --reference POC1605_ref.fa \
        --fastq-dir path/to/barcode10 \
        --out-dir consensus \
        --name barcode10
"""
import argparse
import shlex
import subprocess
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-r", "--reference", required=True, type=Path,
                   help="Reference FASTA to align against and call variants on.")
    p.add_argument("-f", "--fastq-dir", required=True, type=Path,
                   help="Directory of FASTQ reads.")
    p.add_argument("-d", "--out-dir", required=True, type=Path,
                   help="Output directory (created if needed).")
    p.add_argument("-n", "--name", default="sample",
                   help="Basename for output files (default: sample).")
    p.add_argument("--preset", default="map-ont", help="minimap2 -ax preset (default: map-ont).")
    p.add_argument("--minimap2", default="minimap2", help="minimap2 executable.")
    return p


def collect_fastqs(fastq_dir: Path):
    files = sorted(
        list(fastq_dir.glob("*.fastq")) + list(fastq_dir.glob("*.fq")) +
        list(fastq_dir.glob("*.fastq.gz")) + list(fastq_dir.glob("*.fq.gz")))
    if not files:
        raise FileNotFoundError(f"No FASTQ files found in {fastq_dir}")
    return files


def main(argv=None) -> None:
    args = build_parser().parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    bam = args.out_dir / f"{args.name}.sorted.bam"
    vcf = args.out_dir / f"{args.name}.vcf.gz"
    consensus = args.out_dir / f"{args.name}_consensus.fasta"

    fastqs = collect_fastqs(args.fastq_dir)
    fastq_args = " ".join(shlex.quote(str(f)) for f in fastqs)
    ref = shlex.quote(str(args.reference))

    # align + sort
    subprocess.run(
        f'{args.minimap2} -ax {args.preset} {ref} {fastq_args} | samtools sort -o "{bam}"',
        shell=True, check=True, executable="/bin/bash")
    subprocess.run(["samtools", "index", str(bam)], check=True)

    # call variants
    subprocess.run(
        f'bcftools mpileup -Ou -f {ref} "{bam}" | bcftools call -mv -Oz -o "{vcf}"',
        shell=True, check=True, executable="/bin/bash")
    subprocess.run(["bcftools", "index", str(vcf)], check=True)

    # apply consensus
    subprocess.run(
        f'bcftools consensus -f {ref} "{vcf}" > "{consensus}"',
        shell=True, check=True, executable="/bin/bash")

    print(f"Consensus FASTA written to: {consensus}")


if __name__ == "__main__":
    main()
