# HELIOS: HEteroplasmy and Isomer Locator for Organelle Sequences

HELIOS is a lightweight analysis workflow for organelle sequencing data that quantifies heteroplasmic variants and structural isomer proportions directly from raw reads, then packages results into tabular outputs and an interactive report for fast interpretation.

This project provides a Python pipeline to:

- Remap raw FASTQ reads to a circular organelle assembly (BWA-MEM or Minimap2).
- Detect and quantify heteroplasmic SNP/Indel signals from BAM pileup.
- Estimate structural isomer proportions from raw reads using isomer-specific k-mer voting.
- Evaluate plastome assemblies in FASTA-only mode when raw reads are unavailable.
- Generate an interactive HTML report with Plotly.
- Export an isomeric graph in GFA format.

## Inputs

- Circular assembly FASTA (`--fasta`)
- Annotation (`--annotation`, optional: `.gff3` or `.gbk`)
- Raw reads FASTQ (`--fastq`, optional; one or more files; `.gz` supported)

If `--fastq` is omitted, HELIOS runs in assembly-only mode. This mode still writes the HTML report, summary JSON, TSV outputs, and GFA graph, but it does not perform read mapping, heteroplasmy calling, or read-backed isomer quantification. Assembly-only tables and reports use `N/A`-style output for read-derived metrics so absence of measurement is not mistaken for biological zero support.

## Outputs

- `alignments/sample.sorted.bam` and `.bai`
- `results/heteroplasmy.tsv`
- `results/isomer_proportions.tsv`
- `results/isomer_graph.gfa`
- `results/summary.json`
- `report/report.html`

In assembly-only mode, the BAM output is omitted, the report is labeled as candidate plastome isomer evaluation, and read-derived support fields are emitted as not quantified rather than `0`.

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

Assembly-only evaluation:

```bash
helios \
  --fasta input/C2_mq40_final_circular.fasta \
  --output runs/C2 \
  --sample-name C2
```

## External dependencies

For full remapping mode, install:

- `bwa` and `samtools`, or
- `minimap2` and `samtools`

Use `--skip-mapping` to run read-based isomer quantification and report generation without BAM remapping.
Assembly-only mode does not require external aligners because no read mapping is performed.

## Assembly-only validation

HELIOS was validated in assembly-only mode on four plastome assemblies in `input/`:

- `C2` - `Artocarpus integer`
- `D1` - `Durio zibethinus`
- `MELOD1` - `Cucumis melo`
- `STB1` - `Fragaria x ananassa`

All four completed successfully as single circular plastome FASTA inputs. Because each file contained one record and no SSC annotation was supplied, HELIOS generated candidate isomer pairs with its fallback inversion strategy: it keeps the original assembly as Isomer A and creates Isomer B by reverse-complementing the middle 30% of the circular sequence (35%-65% of the contig span). The resulting assembly-only outputs include HTML, TSV, GFA, and summary files without BAM alignments or heteroplasmy calls.
