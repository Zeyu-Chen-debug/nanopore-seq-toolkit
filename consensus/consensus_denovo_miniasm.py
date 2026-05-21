#!/usr/bin/env python3
"""Reference-free assembly with minimap2 all-vs-all overlaps + miniasm.

Merges reads into one FASTQ, computes all-vs-all overlaps, runs miniasm to lay
out contigs, and extracts the contig sequences from the resulting GFA.

Example
-------
    python consensus_denovo_miniasm.py \
        --fastq-dir path/to/barcode21 \
        --out-dir barcode21_consensus_denovo \
        --name barcode21 \
        --minimap2 ./minimap2/minimap2 \
        --miniasm ./miniasm/miniasm
"""
import argparse
import gzip
import shutil
import subprocess
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-f", "--fastq-dir", required=True, type=Path,
                   help="Directory of FASTQ reads.")
    p.add_argument("-d", "--out-dir", required=True, type=Path,
                   help="Output directory (created if needed).")
    p.add_argument("-n", "--name", default="sample",
                   help="Basename for the consensus FASTA (default: sample).")
    p.add_argument("--overlap-mode", default="ava-ont",
                   help="minimap2 all-vs-all preset: ava-ont (Nanopore) or ava-pb (default: ava-ont).")
    p.add_argument("--threads", type=int, default=8, help="minimap2 threads (default: 8).")
    p.add_argument("--minimap2", default="minimap2", help="Path to the minimap2 executable.")
    p.add_argument("--miniasm", default="miniasm", help="Path to the miniasm executable.")
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

    reads_fq = args.out_dir / "reads.fq"
    paf_gz = args.out_dir / "reads.paf.gz"
    gfa = args.out_dir / "reads.gfa"
    consensus = args.out_dir / f"{args.name}_consensus.fasta"

    fastqs = collect_fastqs(args.fastq_dir)
    print(f"Found {len(fastqs)} FASTQ files")

    # merge reads
    with open(reads_fq, "wb") as out_fh:
        for fq in fastqs:
            opener = gzip.open if fq.name.endswith(".gz") else open
            with opener(fq, "rb") as in_fh:
                shutil.copyfileobj(in_fh, out_fh)
    print(f"Merged reads written to: {reads_fq}")

    # all-vs-all overlap
    with open(paf_gz, "wb") as paf_out:
        p1 = subprocess.Popen([args.minimap2, "-x", args.overlap_mode,
                               f"-t{args.threads}", str(reads_fq), str(reads_fq)],
                              stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["gzip", "-1"], stdin=p1.stdout, stdout=paf_out)
        p1.stdout.close()
        r2, r1 = p2.wait(), p1.wait()
        if r1 != 0 or r2 != 0:
            raise RuntimeError(f"overlap step failed (codes {r1}, {r2})")
    print(f"Overlap file written to: {paf_gz}")

    # layout
    with open(gfa, "w") as gfa_out:
        subprocess.run([args.miniasm, "-f", str(reads_fq), str(paf_gz)],
                       check=True, stdout=gfa_out)
    print(f"GFA written to: {gfa}")

    # extract contigs
    contig_num = 0
    with open(consensus, "w") as fasta_out, open(gfa) as gfa_in:
        for line in gfa_in:
            if line.startswith("S\t"):
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                seq = fields[2]
                contig_num += 1
                fasta_out.write(f">contig_{contig_num}\n")
                for i in range(0, len(seq), 80):
                    fasta_out.write(seq[i:i + 80] + "\n")

    if contig_num == 0:
        print("No contigs were found in the GFA output.")
    else:
        print(f"Consensus FASTA written to: {consensus} ({contig_num} contigs)")


if __name__ == "__main__":
    main()
