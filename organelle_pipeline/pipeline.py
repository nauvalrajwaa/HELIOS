"""End-to-end HELIOS pipeline orchestration (FASTQ-required)."""

from __future__ import annotations

from pathlib import Path

from organelle_pipeline.heteroplasmy import call_heteroplasmy_from_bam
from organelle_pipeline.isomer import (
    IsomerPair,
    build_isomer_candidates,
    quantify_isomers_by_remapping,
    quantify_isomers_from_fastq,
    write_isomer_gfa,
)
from organelle_pipeline.mapping import choose_aligner, run_mapping
from organelle_pipeline.models import (
    HETEROPOLASMY_HEADERS,
    ISOMER_HEADERS,
    IsomerQuantResult,
    PipelineConfig,
    PipelineResult,
)
from organelle_pipeline.parsers import (
    estimate_read_length,
    infer_ssc_region,
    read_fasta_records,
)
from organelle_pipeline.repeats import detect_inverted_repeats
from organelle_pipeline.report import write_html_report
from organelle_pipeline.utils import write_json, write_tsv

_MIN_SSC_FRACTION = 0.02  # an SSC must span >=2% of the contig


def _validate_inputs(config: PipelineConfig) -> None:
    if not config.fasta.exists():
        raise FileNotFoundError(f"assembly FASTA not found: {config.fasta}")
    if not config.fastq_files:
        raise ValueError(
            "no FASTQ inputs supplied. HELIOS requires raw reads for every run "
            "(pass one or more files via --fastq; .fastq/.fq and .gz are accepted)"
        )
    missing = [str(p) for p in config.fastq_files if not p.exists()]
    if missing:
        raise FileNotFoundError(f"FASTQ file(s) not found: {', '.join(missing)}")
    if config.annotation is not None and not config.annotation.exists():
        raise FileNotFoundError(f"annotation file not found: {config.annotation}")
    if config.kmer_size < 11:
        raise ValueError("--kmer-size must be >= 11")
    if config.min_depth < 1:
        raise ValueError("--min-depth must be >= 1")
    if config.method not in ("both", "remap", "kmer"):
        raise ValueError(f"unknown quantification method: {config.method}")


def _resolve_ssc_region(config: PipelineConfig, contig_length: int):
    """Return ((start, end), source) or (None, None) per CLI/annotation inputs."""
    if config.ssc_start is not None or config.ssc_end is not None:
        if config.ssc_start is None or config.ssc_end is None:
            raise ValueError("--ssc-start and --ssc-end must be provided together")
        start, end = int(config.ssc_start), int(config.ssc_end)
        if not (0 <= start < end <= contig_length):
            raise ValueError(
                f"--ssc-start/--ssc-end must satisfy 0 <= start < end <= {contig_length}"
            )
        if end - start < max(200, contig_length * _MIN_SSC_FRACTION):
            raise ValueError(
                f"SSC region {start}-{end} spans only {end - start} bp of a "
                f"{contig_length} bp contig — refusing implausible override"
            )
        return (start, end), "user_override"

    if config.annotation is not None:
        region_1based = infer_ssc_region(config.annotation)
        if region_1based:
            start = max(0, region_1based[0] - 1)
            end = min(contig_length, region_1based[1])
            if end - start >= max(200, contig_length * _MIN_SSC_FRACTION):
                return (start, end), "annotation"
    return None, None


def _detect_ssc(sequence: str):
    """Auto-detect IR boundaries with a relaxing length ladder."""
    detection = None
    for min_ir in (10_000, 6_000, 3_000):
        detection = detect_inverted_repeats(sequence, min_ir_len=min_ir)
        if detection.ok:
            return detection
    return detection  # last attempt carries the explanatory warnings


def _build_isomers(
    records: list[tuple[str, str]], config: PipelineConfig
) -> tuple[IsomerPair | None, list[str]]:
    notes: list[str] = []
    sequence = records[0][1]
    region, source = _resolve_ssc_region(config, len(sequence))

    if len(records) == 1 and region is None:
        detection = _detect_ssc(sequence)
        try:
            pair = build_isomer_candidates(records, ssc_region=None, detection=detection)
            notes.append("SSC boundaries auto-detected from self-identical inverted repeats.")
            return pair, notes
        except ValueError as error:
            notes.append(
                f"Isomer quantification skipped: {error} Heteroplasmy-only report generated."
            )
            return None, notes

    ssc_source = source if region else "annotation"
    pair = build_isomer_candidates(records, ssc_region=region, ssc_source=ssc_source)
    return pair, notes


def run_pipeline(config: PipelineConfig) -> PipelineResult:
    """Execute remapping, heteroplasmy calling and isomer quantification."""
    _validate_inputs(config)

    fasta_path: Path = config.fasta
    fastq_paths = [Path(p) for p in config.fastq_files]
    output_dir = Path(config.output_dir)
    alignments_dir = output_dir / "alignments"
    results_dir = output_dir / "results"
    report_dir = output_dir / "report"
    for folder in (output_dir, alignments_dir, results_dir, report_dir):
        folder.mkdir(parents=True, exist_ok=True)

    records = read_fasta_records(fasta_path)
    contig_name, contig_sequence = records[0]
    ref_length = len(contig_sequence)

    aligner = config.aligner
    avg_read_len = estimate_read_length(fastq_paths, limit=config.read_limit_for_stats)
    if aligner == "auto":
        aligner = choose_aligner("auto", avg_read_len)

    bam_a_path = alignments_dir / f"{config.sample_name}.sorted.bam"
    run_mapping(aligner, fasta_path, fastq_paths, bam_a_path, threads=config.threads)

    # ---- isomer references & second alignment ------------------------------
    isomer_pair, assumptions = _build_isomers(records, config)
    bam_b_path: Path | None = None
    b_reference_fasta: Path | None = None
    if isomer_pair is not None and isomer_pair.ssc_start is not None:
        b_reference_fasta = results_dir / "isomer_b_reference.fasta"
        b_reference_fasta.write_text(
            f">{contig_name}\n{isomer_pair.sequence_b}\n", encoding="utf-8"
        )
        bam_b_path = alignments_dir / f"{config.sample_name}.isomerB.sorted.bam"
        run_mapping(aligner, b_reference_fasta, fastq_paths, bam_b_path, threads=config.threads)

    # ---- heteroplasmy -------------------------------------------------------
    heteroplasmy_calls = call_heteroplasmy_from_bam(
        bam_a_path,
        records,
        min_mapq=config.min_mapq,
        min_baseq=config.min_baseq,
        min_depth=config.min_depth,
        min_alt_count=config.min_alt_count,
        min_alt_fraction=config.min_alt_fraction,
    )

    # ---- isomer quantification ----------------------------------------------
    isomer_result: IsomerQuantResult | None = None
    if isomer_pair is not None:
        remap_result = None
        want_remap = config.method in ("both", "remap")
        want_kmer = config.method in ("both", "kmer")

        if want_remap and isomer_pair.ssc_start is not None and bam_b_path is not None:
            remap_result = quantify_isomers_by_remapping(
                bam_a_path,
                bam_b_path,
                contig_name,
                ref_length,
                isomer_pair.ssc_start,
                isomer_pair.ssc_end,
                junction_buffer=config.junction_buffer,
                min_mapq=config.min_mapq,
            )
            if remap_result is None:
                assumptions.append(
                    "Dual-reference remapping produced no informative reads inside the "
                    "SSC core "
                    f"(junction_buffer={config.junction_buffer}); falling back to k-mer voting."
                )

        kmer_result = (
            quantify_isomers_from_fastq(
                fastq_paths,
                isomer_pair,
                kmer_size=config.kmer_size,
                min_hits=config.min_isomer_hits,
                read_limit=config.read_limit_for_isomer,
                assumptions=None,
            )
            if want_kmer
            else None
        )

        isomer_result = _combine_quant_results(
            isomer_pair, remap_result, kmer_result, config.method, assumptions
        )
    else:
        isomer_result = IsomerQuantResult(
            isomer_a_name=contig_name,
            isomer_b_name="NA",
            assigned_a=0,
            assigned_b=0,
            ambiguous=0,
            unassigned=0,
            total_reads_seen=0,
            isomer_a_fraction=0.0,
            isomer_b_fraction=0.0,
            method="not_available",
            assumptions=assumptions,
        )

    # ---- outputs -------------------------------------------------------------
    heteroplasmy_tsv = results_dir / "heteroplasmy.tsv"
    isomer_tsv = results_dir / "isomer_proportions.tsv"
    gfa_graph = results_dir / "isomer_graph.gfa"
    html_report = report_dir / "report.html"
    summary_json = results_dir / "summary.json"

    write_tsv(
        heteroplasmy_tsv,
        [call.to_row() for call in heteroplasmy_calls],
        headers=HETEROPOLASMY_HEADERS,
    )
    write_tsv(isomer_tsv, [isomer_result.to_row()], headers=ISOMER_HEADERS)

    if isomer_pair is not None:
        write_isomer_gfa(gfa_graph, isomer_pair, isomer_result)
    else:
        gfa_graph.write_text("H\tVN:Z:2.0\n", encoding="utf-8")

    write_html_report(html_report, config.sample_name, isomer_result, heteroplasmy_calls)

    result = PipelineResult(
        sample_name=config.sample_name,
        bam_path=bam_a_path,
        heteroplasmy_calls=heteroplasmy_calls,
        isomer_result=isomer_result,
        html_report=html_report,
        gfa_graph=gfa_graph,
        heteroplasmy_tsv=heteroplasmy_tsv,
        isomer_tsv=isomer_tsv,
        summary_json=summary_json,
    )
    write_json(summary_json, result.to_summary())
    return result


def _combine_quant_results(
    isomer_pair: IsomerPair,
    remap_result,
    kmer_result: IsomerQuantResult | None,
    method: str,
    assumptions: list[str],
) -> IsomerQuantResult:
    """Merge estimator outputs into one result record.

    Primary fractions come from dual-reference remapping when available;
    otherwise from k-mer voting.
    """
    base_notes = list(assumptions)
    common = {
        "isomer_a_name": isomer_pair.name_a,
        "isomer_b_name": isomer_pair.name_b,
        "ssc_start": isomer_pair.ssc_start,
        "ssc_end": isomer_pair.ssc_end,
        "ssc_source": isomer_pair.ssc_source,
    }

    if remap_result is not None:
        methods_used = ["dual_reference_remap"]
        if kmer_result is not None:
            methods_used.append("unique_kmer_voting")
        notes = [
            *base_notes,
            "Primary estimate: orientation counting inside the SSC core against both "
            f"isomer references (junction_buffer={remap_result.junction_buffer}bp, "
            f"min_mapq={remap_result.min_mapq}).",
            f"Dual-reference agreement delta: {remap_result.agreement_delta:.4f} "
            "(small values indicate high confidence).",
        ]
        return IsomerQuantResult(
            assigned_a=remap_result.sense_a,
            assigned_b=remap_result.antisense_a,
            ambiguous=kmer_result.ambiguous if kmer_result else 0,
            unassigned=kmer_result.unassigned if kmer_result else 0,
            total_reads_seen=(
                kmer_result.total_reads_seen if kmer_result else remap_result.total_ssc_reads_a
            ),
            isomer_a_fraction=round(remap_result.fraction_a, 6),
            isomer_b_fraction=round(remap_result.fraction_b, 6),
            method="+".join(methods_used),
            assumptions=notes,
            remap_fraction_a=round(remap_result.fraction_a, 6),
            remap_fraction_b=round(remap_result.fraction_b, 6),
            remap_agreement_delta=round(remap_result.agreement_delta, 6),
            remap_total_ssc_reads_a=remap_result.total_ssc_reads_a,
            remap_total_ssc_reads_b=remap_result.total_ssc_reads_b,
            kmer_fraction_a=kmer_result.isomer_a_fraction if kmer_result else None,
            kmer_fraction_b=kmer_result.isomer_b_fraction if kmer_result else None,
            **common,
        )

    if kmer_result is not None:
        merged = kmer_result
        merged.ssc_start = isomer_pair.ssc_start
        merged.ssc_end = isomer_pair.ssc_end
        merged.ssc_source = isomer_pair.ssc_source
        merged.kmer_fraction_a = merged.isomer_a_fraction
        merged.kmer_fraction_b = merged.isomer_b_fraction
        merged.assumptions = base_notes + merged.assumptions[1:]
        return merged

    raise RuntimeError("no isomer quantification estimator produced a result")
