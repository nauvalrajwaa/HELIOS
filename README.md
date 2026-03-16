# HELIOS: HEteroplasmy and Isomer Locator for Organelle Sequences

HELIOS is a lightweight analysis workflow for organelle sequencing data that quantifies heteroplasmic variants and structural isomer proportions directly from raw reads, then packages results into tabular outputs and an interactive report for fast interpretation.

This project provides a Python pipeline to:

- Remap raw FASTQ reads to a circular organelle assembly (BWA-MEM or Minimap2).
- Detect and quantify heteroplasmic SNP/Indel signals from BAM pileup.
- Estimate structural isomer proportions from raw reads using isomer-specific k-mer voting.
- Generate an interactive HTML report with Plotly.
- Export an isomeric graph in GFA format.

## Inputs

- Circular assembly FASTA (`--fasta`)
- Annotation (`--annotation`, optional: `.gff3` or `.gbk`)
- Raw reads FASTQ (`--fastq`, one or more files; `.gz` supported)

## Outputs

- `alignments/sample.sorted.bam` and `.bai`
- `results/heteroplasmy.tsv`
- `results/isomer_proportions.tsv`
- `results/isomer_graph.gfa`
- `results/summary.json`
- `report/report.html`

## Install

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

## Run

```bash
helios \
  --fasta data/organelle.fasta \
  --annotation data/annotation.gff3 \
  --fastq data/reads_R1.fastq.gz data/reads_R2.fastq.gz \
  --output output
```

## External dependencies

For full remapping mode, install:

- `bwa` and `samtools`, or
- `minimap2` and `samtools`

Use `--skip-mapping` to run read-based isomer quantification and report generation without BAM remapping.
