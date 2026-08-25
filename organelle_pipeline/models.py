"""Data models shared across the HELIOS pipeline."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path


@dataclass(slots=True)
class PipelineConfig:
    fasta: Path
    annotation: Path | None
    fastq_files: list[Path]
    output_dir: Path
    sample_name: str
    aligner: str
    threads: int
    min_mapq: int
    min_baseq: int
    min_depth: int
    min_alt_count: int
    min_alt_fraction: float
    kmer_size: int
    min_isomer_hits: int
    read_limit_for_stats: int
    read_limit_for_isomer: int
    method: str = "both"  # both | remap | kmer
    ssc_start: int | None = None
    ssc_end: int | None = None
    junction_buffer: int = 250


@dataclass(slots=True)
class VariantCall:
    contig: str
    position: int
    ref: str
    alt: str
    depth: int
    alt_count: int
    alt_fraction: float
    variant_type: str
    # artifact diagnostics
    strand_ref_fwd: int = 0
    strand_ref_rev: int = 0
    strand_alt_fwd: int = 0
    strand_alt_rev: int = 0
    strand_bias_pvalue: float | None = None
    homopolymer_context: bool = False
    mean_alt_read_position: float | None = None

    def to_row(self) -> dict[str, str]:
        pvalue = f"{self.strand_bias_pvalue:.3g}" if self.strand_bias_pvalue is not None else "NA"
        mean_pos = (
            f"{self.mean_alt_read_position:.3f}"
            if self.mean_alt_read_position is not None
            else "NA"
        )
        return {
            "contig": self.contig,
            "position": str(self.position),
            "ref": self.ref,
            "alt": self.alt,
            "depth": str(self.depth),
            "alt_count": str(self.alt_count),
            "alt_fraction": f"{self.alt_fraction:.6f}",
            "variant_type": self.variant_type,
            "strand_ref_fwd": str(self.strand_ref_fwd),
            "strand_ref_rev": str(self.strand_ref_rev),
            "strand_alt_fwd": str(self.strand_alt_fwd),
            "strand_alt_rev": str(self.strand_alt_rev),
            "strand_bias_p": pvalue,
            "homopolymer": "true" if self.homopolymer_context else "false",
            "mean_alt_read_pos": mean_pos,
        }


HETEROPOLASMY_HEADERS = [
    "contig",
    "position",
    "ref",
    "alt",
    "depth",
    "alt_count",
    "alt_fraction",
    "variant_type",
    "strand_ref_fwd",
    "strand_ref_rev",
    "strand_alt_fwd",
    "strand_alt_rev",
    "strand_bias_p",
    "homopolymer",
    "mean_alt_read_pos",
]


@dataclass(slots=True)
class IsomerQuantResult:
    isomer_a_name: str
    isomer_b_name: str
    assigned_a: int
    assigned_b: int
    ambiguous: int
    unassigned: int
    total_reads_seen: int
    isomer_a_fraction: float
    isomer_b_fraction: float
    method: str
    assumptions: list[str] = field(default_factory=list)
    # dual-reference remapping block (None when unavailable)
    remap_fraction_a: float | None = None
    remap_fraction_b: float | None = None
    remap_agreement_delta: float | None = None
    remap_total_ssc_reads_a: int | None = None
    remap_total_ssc_reads_b: int | None = None
    # k-mer voting block (mirrors fractions when that estimator ran)
    kmer_fraction_a: float | None = None
    kmer_fraction_b: float | None = None
    # SSC provenance
    ssc_start: int | None = None
    ssc_end: int | None = None
    ssc_source: str = "unknown"

    def to_row(self) -> dict[str, str]:
        def fmt(value: float | int | None) -> str:
            if value is None:
                return "NA"
            if isinstance(value, float):
                return f"{value:.6f}"
            return str(value)

        return {
            "isomer_a_name": self.isomer_a_name,
            "isomer_b_name": self.isomer_b_name,
            "assigned_a": str(self.assigned_a),
            "assigned_b": str(self.assigned_b),
            "ambiguous": str(self.ambiguous),
            "unassigned": str(self.unassigned),
            "total_reads_seen": str(self.total_reads_seen),
            "isomer_a_fraction": fmt(self.isomer_a_fraction),
            "isomer_b_fraction": fmt(self.isomer_b_fraction),
            "method": self.method,
            "remap_fraction_a": fmt(self.remap_fraction_a),
            "remap_fraction_b": fmt(self.remap_fraction_b),
            "remap_agreement_delta": fmt(self.remap_agreement_delta),
            "remap_total_ssc_reads_a": fmt(self.remap_total_ssc_reads_a),
            "remap_total_ssc_reads_b": fmt(self.remap_total_ssc_reads_b),
            "kmer_fraction_a": fmt(self.kmer_fraction_a),
            "kmer_fraction_b": fmt(self.kmer_fraction_b),
            "ssc_start": fmt(self.ssc_start),
            "ssc_end": fmt(self.ssc_end),
            "ssc_source": self.ssc_source,
        }


ISOMER_HEADERS = [
    "isomer_a_name",
    "isomer_b_name",
    "assigned_a",
    "assigned_b",
    "ambiguous",
    "unassigned",
    "total_reads_seen",
    "isomer_a_fraction",
    "isomer_b_fraction",
    "method",
    "remap_fraction_a",
    "remap_fraction_b",
    "remap_agreement_delta",
    "remap_total_ssc_reads_a",
    "remap_total_ssc_reads_b",
    "kmer_fraction_a",
    "kmer_fraction_b",
    "ssc_start",
    "ssc_end",
    "ssc_source",
]


@dataclass(slots=True)
class PipelineResult:
    sample_name: str
    bam_path: Path | None
    heteroplasmy_calls: list[VariantCall]
    isomer_result: IsomerQuantResult
    html_report: Path
    gfa_graph: Path
    heteroplasmy_tsv: Path
    isomer_tsv: Path
    summary_json: Path

    def to_summary(self) -> dict:
        return {
            "sample_name": self.sample_name,
            "bam_path": str(self.bam_path) if self.bam_path else None,
            "heteroplasmy_call_count": len(self.heteroplasmy_calls),
            "isomer": asdict(self.isomer_result),
            "output_files": {
                "bam": str(self.bam_path) if self.bam_path else None,
                "report_html": str(self.html_report),
                "isomer_graph_gfa": str(self.gfa_graph),
                "heteroplasmy_tsv": str(self.heteroplasmy_tsv),
                "isomer_proportions_tsv": str(self.isomer_tsv),
                "summary_json": str(self.summary_json),
            },
        }
