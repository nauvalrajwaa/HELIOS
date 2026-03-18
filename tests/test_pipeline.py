import json
from pathlib import Path

from organelle_pipeline.models import PipelineConfig
from organelle_pipeline.pipeline import run_pipeline


def test_run_pipeline_supports_assembly_only_mode(tmp_path: Path):
    fasta_path = tmp_path / "assembly.fasta"
    fasta_path.write_text(
        ">plastome\n"
        "ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"
        "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n"
    )

    output_dir = tmp_path / "out"
    config = PipelineConfig(
        fasta=fasta_path,
        annotation=None,
        fastq_files=[],
        output_dir=output_dir,
        sample_name="assembly-only",
        aligner="auto",
        threads=1,
        min_mapq=20,
        min_baseq=20,
        min_depth=200,
        min_alt_count=10,
        min_alt_fraction=0.01,
        kmer_size=31,
        min_isomer_hits=2,
        skip_mapping=True,
        read_limit_for_stats=100,
        read_limit_for_isomer=100,
    )

    result = run_pipeline(config)

    assert result.bam_path is None
    assert result.heteroplasmy_calls == []
    assert result.isomer_result.method == "assembly_only_candidates"
    assert result.isomer_result.total_reads_seen == 0
    assert result.heteroplasmy_tsv.exists()
    assert result.isomer_tsv.exists()
    assert result.gfa_graph.exists()
    assert result.html_report.exists()
    assert result.summary_json.exists()

    summary = json.loads(result.summary_json.read_text())
    assert summary["bam_path"] is None
    assert summary["heteroplasmy_call_count"] == 0
    assert summary["isomer"]["method"] == "assembly_only_candidates"

    report_html = result.html_report.read_text()
    assert "Assembly-only evaluation" in report_html
    assert "candidate plastome isomers without read-backed proportions" in report_html
    assert "Reads Processed</div><div class=\"value\">N/A" in report_html
    assert "Heteroplasmy calling is disabled in assembly-only mode." in report_html

    isomer_tsv = result.isomer_tsv.read_text()
    assert "NA" in isomer_tsv
    assert "assembly_only_candidates" in isomer_tsv

    gfa_graph = result.gfa_graph.read_text()
    assert "CL:Z:assembly_only_candidates" in gfa_graph
    assert "RC:i:" not in gfa_graph
    assert "PR:f:" not in gfa_graph
