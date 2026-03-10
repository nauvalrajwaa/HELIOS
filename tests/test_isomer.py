from pathlib import Path

from organelle_pipeline.isomer import build_isomer_candidates
from organelle_pipeline.isomer import quantify_isomers_from_fastq


def test_build_isomer_candidates_from_two_records():
    records = [("isoA", "ACGTACGTACGT"), ("isoB", "ACGTTCGTACGA")]
    a, b, assumptions = build_isomer_candidates(records, None)
    assert a[0] == "isoA"
    assert b[0] == "isoB"
    assert assumptions


def test_quantify_isomers_from_fastq(tmp_path: Path):
    fastq_path = tmp_path / "reads.fastq"
    fastq_path.write_text(
        "@r1\nAAAAAACCCCC\n+\nFFFFFFFFFFF\n"
        "@r2\nAAAAAACCCCC\n+\nFFFFFFFFFFF\n"
        "@r3\nTTTTTTGGGGG\n+\nFFFFFFFFFFF\n"
        "@r4\nTTTTTTGGGGG\n+\nFFFFFFFFFFF\n",
        encoding="utf-8",
    )

    result = quantify_isomers_from_fastq(
        fastq_files=[fastq_path],
        isomer_a=("A", "AAAAAACCCCCAAAAAACCCCC"),
        isomer_b=("B", "GGGGGTTTTTGGGGGTTTTT"),
        kmer_size=5,
        min_hits=1,
        read_limit=100,
        assumptions=["test"],
    )

    assert result.assigned_a == 2
    assert result.assigned_b == 2
    assert result.total_reads_seen == 4
    assert result.isomer_a_fraction == 0.5
    assert result.isomer_b_fraction == 0.5
