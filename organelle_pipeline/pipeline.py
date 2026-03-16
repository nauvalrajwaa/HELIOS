from pathlib import Path

from organelle_pipeline.heteroplasmy import call_heteroplasmy_from_bam
from organelle_pipeline.isomer import build_isomer_candidates
from organelle_pipeline.isomer import quantify_isomers_from_fastq
from organelle_pipeline.isomer import write_isomer_gfa
from organelle_pipeline.mapping import choose_aligner
from organelle_pipeline.mapping import run_mapping
from organelle_pipeline.models import PipelineConfig
from organelle_pipeline.models import PipelineResult
from organelle_pipeline.parsers import estimate_read_length
from organelle_pipeline.parsers import infer_ssc_region
from organelle_pipeline.parsers import read_fasta_records
from organelle_pipeline.report import write_html_report
from organelle_pipeline.utils import write_json
from organelle_pipeline.utils import write_tsv


def run_pipeline(config: PipelineConfig) -> PipelineResult:
    _validate_inputs(config)

    fasta_records = read_fasta_records(config.fasta)
    average_read_length = estimate_read_length(config.fastq_files, limit=config.read_limit_for_stats)
    aligner = choose_aligner(config.aligner, average_read_length)

    alignments_dir = config.output_dir / "alignments"
    results_dir = config.output_dir / "results"
    report_dir = config.output_dir / "report"
    alignments_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)

    bam_path: Path | None = None
    if not config.skip_mapping:
        bam_path = alignments_dir / f"{config.sample_name}.sorted.bam"
        run_mapping(
            aligner=aligner,
            fasta_path=config.fasta,
            fastq_files=config.fastq_files,
            output_bam=bam_path,
            threads=config.threads,
        )

    if bam_path:
        heteroplasmy_calls = call_heteroplasmy_from_bam(
            bam_path=bam_path,
            fasta_records=fasta_records,
            min_mapq=config.min_mapq,
            min_baseq=config.min_baseq,
            min_depth=config.min_depth,
            min_alt_count=config.min_alt_count,
            min_alt_fraction=config.min_alt_fraction,
        )
    else:
        heteroplasmy_calls = []

    ssc_region = infer_ssc_region(config.annotation)
    isomer_a, isomer_b, assumptions = build_isomer_candidates(fasta_records, ssc_region)
    isomer_result = quantify_isomers_from_fastq(
        fastq_files=config.fastq_files,
        isomer_a=isomer_a,
        isomer_b=isomer_b,
        kmer_size=config.kmer_size,
        min_hits=config.min_isomer_hits,
        read_limit=config.read_limit_for_isomer,
        assumptions=assumptions,
    )

    heteroplasmy_tsv = results_dir / "heteroplasmy.tsv"
    isomer_tsv = results_dir / "isomer_proportions.tsv"
    gfa_graph = results_dir / "isomer_graph.gfa"
    report_html = report_dir / "report.html"
    summary_json = results_dir / "summary.json"

    write_tsv(
        heteroplasmy_tsv,
        rows=[call.to_row() for call in heteroplasmy_calls],
        headers=["contig", "position", "ref", "alt", "depth", "alt_count", "alt_fraction", "variant_type"],
    )
    write_tsv(
        isomer_tsv,
        rows=[isomer_result.to_row()],
        headers=[
            "isomer_a",
            "isomer_b",
            "assigned_a",
            "assigned_b",
            "ambiguous",
            "unassigned",
            "total_reads_seen",
            "isomer_a_fraction",
            "isomer_b_fraction",
            "method",
        ],
    )
    write_isomer_gfa(gfa_graph, isomer_result)
    write_html_report(report_html, config.sample_name, isomer_result, heteroplasmy_calls)

    result = PipelineResult(
        sample_name=config.sample_name,
        bam_path=bam_path,
        heteroplasmy_calls=heteroplasmy_calls,
        isomer_result=isomer_result,
        html_report=report_html,
        gfa_graph=gfa_graph,
        heteroplasmy_tsv=heteroplasmy_tsv,
        isomer_tsv=isomer_tsv,
        summary_json=summary_json,
    )
    write_json(summary_json, result.to_summary())
    return result


def _validate_inputs(config: PipelineConfig) -> None:
    if not config.fasta.exists():
        raise FileNotFoundError(f"Missing FASTA: {config.fasta}")
    for fastq_file in config.fastq_files:
        if not fastq_file.exists():
            raise FileNotFoundError(f"Missing FASTQ: {fastq_file}")
    if config.annotation and not config.annotation.exists():
        raise FileNotFoundError(f"Missing annotation: {config.annotation}")
    if config.kmer_size < 11:
        raise ValueError("kmer_size must be >= 11 for stable isomer voting")
    if config.min_depth < 1:
        raise ValueError("min_depth must be >= 1")


if __name__ == "__main__":
    from organelle_pipeline.__main__ import main

    raise SystemExit(main())
