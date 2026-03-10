from dataclasses import asdict
from dataclasses import dataclass
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
    skip_mapping: bool
    read_limit_for_stats: int
    read_limit_for_isomer: int


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

    def to_row(self) -> dict[str, str | int | float]:
        return {
            "contig": self.contig,
            "position": self.position,
            "ref": self.ref,
            "alt": self.alt,
            "depth": self.depth,
            "alt_count": self.alt_count,
            "alt_fraction": round(self.alt_fraction, 6),
            "variant_type": self.variant_type,
        }


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
    assumptions: list[str]

    def to_row(self) -> dict[str, str | int | float]:
        return {
            "isomer_a": self.isomer_a_name,
            "isomer_b": self.isomer_b_name,
            "assigned_a": self.assigned_a,
            "assigned_b": self.assigned_b,
            "ambiguous": self.ambiguous,
            "unassigned": self.unassigned,
            "total_reads_seen": self.total_reads_seen,
            "isomer_a_fraction": round(self.isomer_a_fraction, 6),
            "isomer_b_fraction": round(self.isomer_b_fraction, 6),
            "method": self.method,
        }


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

    def to_summary(self) -> dict[str, object]:
        return {
            "sample_name": self.sample_name,
            "bam_path": str(self.bam_path) if self.bam_path else None,
            "heteroplasmy_call_count": len(self.heteroplasmy_calls),
            "isomer": asdict(self.isomer_result),
            "output_files": {
                "html_report": str(self.html_report),
                "gfa_graph": str(self.gfa_graph),
                "heteroplasmy_tsv": str(self.heteroplasmy_tsv),
                "isomer_tsv": str(self.isomer_tsv),
                "summary_json": str(self.summary_json),
            },
        }
