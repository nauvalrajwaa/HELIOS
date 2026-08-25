# HELIOS: HEteroplasmy and Isomer Locator for Organelle Sequences

HELIOS is a lightweight analysis workflow for organelle sequencing data that quantifies heteroplasmic variants and structural isomer proportions directly from raw reads, then packages results into tabular outputs and an interactive report for fast interpretation.

This project provides a Python pipeline to:

- Remap raw FASTQ reads to a circular organelle assembly (BWA-MEM or Minimap2).
- Detect and quantify heteroplasmic SNP/Indel signals from BAM pileup with strand-bias statistics, VCF-style indel alleles, homopolymer-context and read-position artifact flags.
- Determine SSC boundaries from the assembly itself via inverted-repeat auto-detection (no annotation required), or accept curated `--ssc-start/--ssc-end` / GFF3/GBK coordinates.
- Estimate isomer proportions by dual-reference remapping (orientation of first-in-pair reads inside the SSC core against both Isomer A and Isomer B references), cross-checked by unique k-mer voting.
- Generate an interactive HTML report with Plotly.
- Export an isomeric graph in GFA format carrying real segment sequences.

## Inputs

- Circular assembly FASTA (`--fasta`, required)
- Raw reads FASTQ (`--fastq`, **required**: one or more files; `.gz` supported)
- Annotation (`--annotation`, optional: `.gff3` or `.gbk`) — used as an SSC hint when provided
- `--ssc-start/--ssc-end` — optional explicit SSC coordinates (0-based half-open); override everything else
- `--method both|remap|kmer` — which estimator(s) to run (default `both`)

Raw reads are mandatory: every read-derived metric in HELIOS is a measurement, and the pipeline refuses to fabricate candidate numbers without them.

## How isomers are resolved

1. **Two-record FASTA** → records are taken directly as Isomer A/B.
2. **Single record + explicit/annotated SSC** → Isomer B = reverse complement of the SSC arc.
3. **Single record, no hint** → HELIOS detects the two inverted-repeat copies (IRA/IRB) itself — k-mer seeding, maximal extension (origin-aware), clustering of partial repeat hits, then plausibility-scored pair selection — and derives the SSC arc between them.

For D1 (`Durio zibethinus`) no credible IR pair exists (reduced/diverged IR), so manual coordinates must be supplied; see `input/ssc_annotations.tsv`.

## Outputs

- `alignments/sample.sorted.bam` + `.bai` (vs Isomer A reference)
- `alignments/sample.isomerB.sorted.bam` + `.bai` (vs Isomer B reference)
- `results/isomer_b_reference.fasta`
- `results/heteroplasmy.tsv` (now with per-call strand counts, Fisher strand-bias p-value, homopolymer flag, mean alt-read position)
- `results/isomer_proportions.tsv` (remap + k-mer estimates, agreement delta, SSC provenance)
- `results/isomer_graph.gfa` (real sequences, assigned-read tags)
- `results/summary.json`
- `report/report.html`

## Install

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e '.[dev]'
```

## Run

```bash
helios \
  --fasta input/C2_mq40_final_circular.fasta \
  --fastq reads_C2_R1.fastq.gz reads_C2_R2.fastq.gz \
  --output runs/plastome/C2 \
  --sample-name C2
```

Useful extras: `--method remap`, `--junction-buffer 250` (bases trimmed from each SSC edge before orientation counting), `--ssc-start/--ssc-end` overrides.

Curate SSC coordinates for a folder of assemblies ahead of time:

```bash
python scripts/curate_ssc_annotations.py input -o input/ssc_annotations.tsv
```

## External dependencies

- `bwa` + `samtools`, or `minimap2` + `samtools` (auto-selected by read length: ≥800 bp average → minimap2 ONT preset, otherwise BWA-MEM)

## Tests

```bash
pytest
```

18 tests cover the statistics helpers, indel allele construction, IR detection on synthetic genomes, k-mer voting (including reverse-complement reads), and an end-to-end run that simulates paired reads with a known 30% Isomer B mixture and asserts the pipeline recovers it.

## Validation status

- Synthetic end-to-end mixture recovery: PASS (30% ± tolerance)
- Assembly self-detection validated against known plastome architecture:
  - C2 (*Artocarpus integer*, 158,855 bp): IR 24,470/24,436 · SSC 19,397
  - MELOD1 (*Cucumis melo*, 155,999 bp): IR 25,794/25,795 · SSC 18,085 (consistent with NCBI OR643681.1: IR 24,114 · SSC 20,963)
  - STB1 (*Fragaria × ananassa*, 155,568 bp): IR 25,935/25,936 · SSC 18,145
  - D1 (*Durio zibethinus*, 142,917 bp): no detectable perfect IR pair — supply SSC manually
- Read-backed quantification on real FASTQ: pending data transfer from HPC

## Project structure

```
organelle_pipeline/
  __main__.py      # helios CLI
  pipeline.py      # orchestration: SSC resolution -> mapping -> calling -> outputs
  repeats.py       # inverted-repeat detection & SSC boundary derivation
  isomer.py        # isomer pair construction, dual-reference remapping, k-mer voting
  heteroplasmy.py  # pileup caller with artifact diagnostics (Fisher strand-bias, etc.)
  mapping.py       # bwa/minimap2 + samtools wrappers with index caching
  parsers.py       # FASTA/FASTQ/GFF/GenBank readers
  models.py        # dataclasses & TSV schemas
  report.py        # interactive HTML report
  utils.py         # shared helpers
scripts/
  curate_ssc_annotations.py  # batch SSC curation for a folder of assemblies
input/              # circular assemblies + curated ssc_annotations.tsv
tests/              # pytest suite incl. end-to-end synthetic mixture recovery
```

## Development

```bash
pip install -e '.[dev]'
ruff check --fix .   # lint (config in pyproject.toml)
ruff format .
pytest
```
