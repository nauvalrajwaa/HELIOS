from collections import Counter

from organelle_pipeline.heteroplasmy import _calls_from_counts


def test_calls_from_counts_returns_snp_and_indel():
    calls = _calls_from_counts(
        contig="chrO",
        position=123,
        ref_base="A",
        base_counter=Counter({"A": 70, "G": 20}),
        indel_counter=Counter({"DEL": 10}),
        depth=100,
        min_depth=30,
        min_alt_count=5,
        min_alt_fraction=0.05,
    )

    assert len(calls) == 2
    snp = [call for call in calls if call.variant_type == "SNP"][0]
    indel = [call for call in calls if call.variant_type == "INDEL"][0]
    assert snp.alt == "G"
    assert round(snp.alt_fraction, 4) == 0.2
    assert indel.alt == "DEL"
    assert round(indel.alt_fraction, 4) == 0.1


def test_calls_from_counts_filters_low_support():
    calls = _calls_from_counts(
        contig="chrO",
        position=10,
        ref_base="C",
        base_counter=Counter({"C": 95, "T": 3}),
        indel_counter=Counter(),
        depth=98,
        min_depth=30,
        min_alt_count=5,
        min_alt_fraction=0.05,
    )
    assert calls == []
