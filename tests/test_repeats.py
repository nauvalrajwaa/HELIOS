"""Unit tests for inverted-repeat / SSC boundary detection."""

import random

from organelle_pipeline.repeats import detect_inverted_repeats, invert_segment
from organelle_pipeline.utils import reverse_complement


def _random_seq(rng: random.Random, n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_detects_synthetic_ir_pair_and_derives_ssc():
    rng = random.Random(42)
    ir = _random_seq(rng, 8_000)
    ssc = _random_seq(rng, 6_000)
    lsc = _random_seq(rng, 20_000)
    genome = lsc + ir + ssc + reverse_complement(ir)
    n = len(genome)

    detection = detect_inverted_repeats(genome)
    assert detection.ok
    ssc_len = (detection.ssc_end - detection.ssc_start) % n
    assert abs(ssc_len - 6_000) <= 4
    assert abs(detection.ir_a.length - detection.ir_b.length) <= 8
    total = (
        detection.ir_a.length
        + detection.ir_b.length
        + ssc_len
        + ((detection.lsc_end - detection.lsc_start) % n)
    )
    assert total == n


def test_handles_ir_copy_spanning_the_assembly_origin():
    rng = random.Random(7)
    head = _random_seq(rng, 8_000)
    ir = _random_seq(rng, 6_000)
    ssc = _random_seq(rng, 4_000)
    tail = _random_seq(rng, 7_000)
    # circle order: head | irA | ssc | tail | irB -- both single-copy arcs > 0,
    # unlike a naive ir+ssc+tail+rc(ir) chain whose copies end up back-to-back
    # through the origin (a palindromic block, not a valid plastome layout).
    base = head + ir + ssc + tail + reverse_complement(ir)
    # rotate so the linear FASTA cut lands inside the second IR copy
    cut = len(head) + len(ir) + len(ssc) + len(tail) + 3_000
    rotated = base[cut:] + base[:cut]

    detection = detect_inverted_repeats(rotated, min_ir_len=4_000)
    assert detection.ok


def test_reports_warning_when_no_ir_pair_exists():
    rng = random.Random(99)
    genome = _random_seq(rng, 30_000)

    detection = detect_inverted_repeats(genome)
    assert not detection.ok
    assert detection.warnings


def test_double_inversion_restores_original_for_non_wrapping_ssc():
    rng = random.Random(5)
    ir = _random_seq(rng, 5_000)
    ssc = _random_seq(rng, 3_000)
    genome = _random_seq(rng, 10_000) + ir + ssc + reverse_complement(ir)

    detection = detect_inverted_repeats(genome, min_ir_len=3_000)
    assert detection.ok
    assert detection.ssc_start < detection.ssc_end  # non-wrapping in this layout

    flipped = invert_segment(genome, detection.ssc_start, detection.ssc_end)
    restored = invert_segment(flipped, detection.ssc_start, detection.ssc_end)
    assert restored == genome
