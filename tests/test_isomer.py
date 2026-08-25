"""Unit tests for isomer construction and both quantification estimators."""

import random

import pytest

from organelle_pipeline.isomer import (
    IsomerPair,
    _build_unique_kmers,
    _ssc_core_pieces,
    build_isomer_candidates,
    quantify_isomers_from_fastq,
)
from organelle_pipeline.repeats import IrDetection, RepeatRegion
from organelle_pipeline.utils import reverse_complement


def test_build_isomer_candidates_from_two_records():
    records = [("isoA", "ACGTACGTACGT"), ("isoB", "ACGTTCGTACGA")]
    pair = build_isomer_candidates(records)
    assert pair.name_a == "isoA"
    assert pair.sequence_a == "ACGTACGTACGT"
    assert pair.name_b == "isoB"
    assert pair.sequence_b == "ACGTTCGTACGA"
    assert pair.ssc_start is None  # no SSC frame for the two-record shortcut
    assert pair.assumptions


def test_build_isomer_candidates_with_ssc_override():
    rng = random.Random(11)
    sequence = "".join(rng.choice("ACGT") for _ in range(800))
    records = [("ctg", sequence)]
    pair = build_isomer_candidates(records, ssc_region=(250, 550), ssc_source="user_override")
    assert pair.ssc_start == 250 and pair.ssc_end == 550
    assert pair.ssc_source == "user_override"
    assert pair.contig_name == "ctg"
    expected = sequence[:250] + reverse_complement(sequence[250:550]) + sequence[550:]
    assert pair.sequence_b == expected


def test_build_isomer_candidates_requires_ssc_for_single_record():
    with pytest.raises(ValueError):
        build_isomer_candidates([("solo", "ACGT" * 50)], detection=None)


def test_auto_detection_via_injection():
    rng = random.Random(12)
    ir = "".join(rng.choice("ACGT") for _ in range(400))
    ssc = "".join(rng.choice("ACGT") for _ in range(300))
    genome = "".join(rng.choice("ACGT") for _ in range(900)) + ir + ssc + reverse_complement(ir)
    detection = IrDetection(
        ir_a=RepeatRegion("IRA", 1300, 1700, 400),
        ir_b=RepeatRegion("IRB", 900, 1300, 400),
        ssc_start=900 + 0,  # arc between copies
        ssc_end=1300,
        lsc_start=None,
        lsc_end=None,
    )
    detection.ssc_start = 900
    detection.ssc_end = 1300
    pair = build_isomer_candidates([("ctg", genome)], detection=detection)
    assert pair.ssc_source == "auto_detected"
    assert pair.ssc_start == 900 and pair.ssc_end == 1300
    expected = genome[:900] + reverse_complement(genome[900:1300]) + genome[1300:]
    assert pair.sequence_b == expected


def test_unique_kmers_symmetric_difference():
    kmers_a, kmers_b = _build_unique_kmers("AAAACCCC", "AAAAGGGG", kmer_size=4)
    assert "AAAA" not in kmers_a and "AAAA" not in kmers_b  # shared k-mer removed
    assert {"AAAC", "AACC", "ACCC", "CCCC"} <= kmers_a
    assert {"AAAG", "AAGG", "AGGG", "GGGG"} <= kmers_b


def test_ssc_core_pieces_non_wrapping():
    pieces = _ssc_core_pieces(1000, 5000, 10_000, 100)
    assert pieces == [(1100, 4900)]


def test_ssc_core_pieces_wrapping():
    pieces = _ssc_core_pieces(9_000, 1_000, 10_000, 100)
    assert pieces == [(9_100, 10_000), (0, 900)]


def test_quantify_isomers_from_fastq_kmer_voting(tmp_path):
    rng = random.Random(3)
    seq_a = "".join(rng.choice("ACGT") for _ in range(80))
    seq_b = "".join(rng.choice("ACGT") for _ in range(80))
    isomer_pair = IsomerPair(
        name_a="A",
        sequence_a=seq_a,
        name_b="B",
        sequence_b=seq_b,
    )

    reads = [
        seq_a[:30],  # A molecule, forward frame
        seq_a[30:60],  # A molecule, forward frame
        seq_b[:30],  # B molecule, forward frame
        reverse_complement(seq_b[30:60]),  # B molecule sequenced from the other strand
        "ACGTACGTACGTACGTACGT",  # unrelated -> unassigned
    ]
    fastq = tmp_path / "reads.fastq"
    lines = []
    for i, read in enumerate(reads):
        lines += [f"@r{i}", read, "+", "I" * len(read)]
    fastq.write_text("\n".join(lines) + "\n")

    result = quantify_isomers_from_fastq(
        [fastq],
        isomer_pair,
        kmer_size=21,
        min_hits=3,
        read_limit=100,
        assumptions=None,
    )
    assert result.assigned_a == 2
    assert result.assigned_b == 2
    assert result.unassigned == 1
    assert result.ambiguous == 0
    assert result.total_reads_seen == 5
    assert result.isomer_a_fraction == pytest.approx(0.5)
    assert result.isomer_b_fraction == pytest.approx(0.5)
