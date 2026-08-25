"""Unit tests for the hardened heteroplasmy caller."""

from collections import Counter

from organelle_pipeline.heteroplasmy import (
    _adjacent_homopolymer,
    _calls_from_counts,
    fisher_exact_2x2,
)


def test_fisher_exact_balanced_table_has_no_strand_bias():
    # alt 5F/5R vs ref 5F/5R -> perfectly balanced
    pvalue = fisher_exact_2x2(5, 5, 5, 5)
    assert pvalue > 0.99


def test_fisher_exact_one_sided_alt_support_flags_bias():
    # all alt reads on forward strand only -> strong bias
    pvalue = fisher_exact_2x2(10, 0, 5, 5)
    assert pvalue < 0.05


def test_fisher_exact_extremes():
    assert fisher_exact_2x2(0, 0, 0, 0) == 1.0
    # overwhelming signal
    assert fisher_exact_2x2(50, 0, 50, 50) < 1e-6


def test_adjacent_homopolymer_detection():
    sequence = "ACGT" + "A" * 10 + "CGTA"
    assert _adjacent_homopolymer(sequence, 6) is True
    plain = "ACGTACGTACGTACGT"
    assert _adjacent_homopolymer(plain, 5) is False


def _stats_bundle(snp_counts=None, snp_strand=None, indel_counts=None, indel_strand=None):
    return (
        "A",  # ref_base
        100,  # depth
        snp_counts or Counter(),
        (40, 40),  # strand_ref fwd/rev
        snp_strand or {},
        {},
        indel_counts or Counter(),
        indel_strand or {},
    )


def test_snp_call_carries_strand_statistics():
    calls = _calls_from_counts(
        "ctg",
        123,
        "A" * 200,
        *_stats_bundle(
            snp_counts=Counter({"G": 20}),
            snp_strand={"G": [18, 2]},
        ),
        min_depth=30,
        min_alt_count=5,
        min_alt_fraction=0.05,
    )
    assert len(calls) == 1
    call = calls[0]
    assert (call.contig, call.position, call.ref, call.alt) == ("ctg", 123, "A", "G")
    assert call.variant_type == "SNP"
    assert call.strand_alt_fwd == 18 and call.strand_alt_rev == 2
    assert call.strand_ref_fwd == 40 and call.strand_ref_rev == 40
    assert call.strand_bias_pvalue is not None and call.strand_bias_pvalue < 0.05


def test_insertion_gets_vcf_style_allele_with_sequence():
    reference = "A" * 130
    calls = _calls_from_counts(
        "ctg",
        123,
        reference,
        *_stats_bundle(
            indel_counts=Counter({("INS", 3, "TTT"): 12}),
            indel_strand={("INS", 3, "TTT"): [7, 5]},
        ),
        min_depth=30,
        min_alt_count=5,
        min_alt_fraction=0.05,
    )
    assert len(calls) == 1
    call = calls[0]
    assert call.variant_type == "INS"
    assert call.ref == "A"
    assert call.alt == "ATTT"


def test_deletion_carries_deleted_reference_sequence():
    reference = "A" * 130
    deleted_seq = reference[124 : 124 + 4]  # 4 bp after anchor position 123
    assert deleted_seq == "AAAA"
    calls = _calls_from_counts(
        "ctg",
        123,
        reference,
        *_stats_bundle(
            indel_counts=Counter({("DEL", 4, deleted_seq): 15}),
            indel_strand={("DEL", 4, deleted_seq): [9, 6]},
        ),
        min_depth=30,
        min_alt_count=5,
        min_alt_fraction=0.05,
    )
    assert len(calls) == 1
    call = calls[0]
    assert call.variant_type == "DEL"
    assert call.ref == "AAAAA"  # anchor + 4 deleted bases
    assert call.alt == "A"


def test_low_depth_and_weak_alleles_are_suppressed():
    calls = _calls_from_counts(
        "ctg",
        123,
        "A" * 200,
        *_stats_bundle(
            snp_counts=Counter({"T": 3}),
            snp_strand={"T": [2, 1]},
        ),
        min_depth=30,
        min_alt_count=5,
        min_alt_fraction=0.05,
    )
    assert calls == []
