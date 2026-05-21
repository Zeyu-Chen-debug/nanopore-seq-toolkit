# Nanopore Sequence Analysis Toolkit

A collection of command-line Python tools for downstream analysis of Oxford
Nanopore long-read sequencing data. They were written to validate large
de-novo-assembled DNA constructs and synthetic genomes: aligning reads to a
reference, assessing read quality and error rates, plotting coverage, locating
mismatches, calling methylation, and building consensus sequences.

Every tool is a self-contained script with an `argparse` command-line
interface — run any of them with `-h` to see its options.

## Installation

```bash
git clone https://github.com/<your-username>/nanopore-seq-toolkit.git
cd nanopore-seq-toolkit
pip install -r requirements.txt
```

### External tools

These scripts orchestrate standard bioinformatics binaries, which must be
installed separately and available on your `PATH` (or passed explicitly, e.g.
`--minimap2 /path/to/minimap2`):

| Tool | Used by |
|------|---------|
| [minimap2](https://github.com/lh3/minimap2) | alignment, coverage, mismatch, consensus |
| [miniasm](https://github.com/lh3/miniasm) | `consensus/consensus_denovo_miniasm.py` |
| [samtools](https://www.htslib.org/) | consensus tools |
| [bcftools](https://www.htslib.org/) | consensus tools |
| [modkit](https://github.com/nanoporetech/modkit) | methylation tools |

## Tools

### `alignment/`
- **`generate_paf.py`** — align reads to a reference with minimap2 and write a
  PAF; optionally export an Excel summary with per-alignment identity, span,
  mismatch count and error rate.
- **`percent_aligned.py`** — report the percentage of reads that aligned
  (total FASTQ reads vs. unique aligned reads in the PAF).

### `read_quality/`
- **`good_reads_summary.py`** — for one or more samples, compute the fraction
  of "good" reads (those covering ≥ a threshold of the target) and produce a
  3-panel figure (good-read %, read-length distribution, error rate) plus an
  Excel table.
- **`average_error_rate.py`** — print the mean error rate per group from an
  alignment summary spreadsheet.
- **`error_rate_boxplot.py`** — box-and-strip plot of per-read error rate
  across labelled alignment groups.
- **`read_length_density.py`** — kernel-density plot of read-length
  distributions across labelled groups (log axis).
- **`good_read_vs_length.py`** — scatter of good-read percentage vs. template
  length, with a linear trendline.
- **`good_read_barplot.py`** — grouped bar chart of a per-replicate metric,
  faceted by template.
- **`summary_boxplots.py`** — grid of box-and-strip plots summarising QC
  metrics by treatment and time point.

### `coverage/`
- **`coverage_plot.py`** — per-base coverage depth across a reference for one
  sample.
- **`insert_coverage_compare.py`** — overlay coverage from several samples,
  with random subsampling for fair comparison at equal depth.
- **`alignment_spans.py`** — draw each read's alignment span along the
  reference (useful for spotting linearisation/cut sites).

### `mismatch/`
- **`mismatch_positions.py`** — stream reads through minimap2 with the `--cs`
  tag, parse substitutions from passing alignments, and plot a per-position
  mismatch count.

### `methylation/`
- **`methylation_pileup.py`** — run `modkit pileup` for a chosen modified base
  (e.g. 5mC, 6mA) and export filtered methylation calls to Excel.
- **`filter_methylation.py`** — filter an existing methylation Excel table by a
  minimum percent-modified value.

### `consensus/`
- **`consensus_reference.py`** — reference-guided consensus
  (minimap2 → samtools → bcftools).
- **`consensus_denovo_polish.py`** — reference-free consensus by seeding on the
  longest read and iteratively polishing.
- **`consensus_denovo_miniasm.py`** — reference-free assembly via minimap2
  all-vs-all overlaps + miniasm.

## Example

```bash
# 1. Align and summarise
python alignment/generate_paf.py \
    --reference ref.fa --fastq-dir reads/barcode01 \
    --output-paf barcode01.paf --to-xlsx

# 2. Compare read quality across conditions
python read_quality/good_reads_summary.py \
    --reference ref.fa \
    --sample "FK=reads/barcode05" \
    --sample "FKS=reads/barcode06" \
    --out-prefix results/comparison

# 3. Build a reference-guided consensus
python consensus/consensus_reference.py \
    --reference ref.fa --fastq-dir reads/barcode10 \
    --out-dir consensus --name barcode10
```

## License

Released under the MIT License — see [LICENSE](LICENSE).
