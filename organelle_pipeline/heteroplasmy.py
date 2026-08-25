"""Heteroplasmy calling from BAM alignments.

Pure-Python pileup walker built on pysam. Emits SNPs plus VCF-style anchored
indel alleles with real inserted/deleted sequences. Every call carries
strand-level support counts, an exact two-sided Fisher test for strand bias,
a homopolymer-context flag and the mean relative alt-read position — the
artifact indicators that matter most on organellar genomes.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from math import exp, lgamma
from pathlib import Path

from organelle_pipeline.models import VariantCall

_HOMOPOLYMER_RUN = 6  # run length considered homopolymer-prone
_HOMOPOLYMER_WINDOW = 10


# ---------------------------------------------------------------------------
# statistics helpers
# ---------------------------------------------------------------------------


def _log_binom(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)


def fisher_exact_2x2(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher exact p-value for table [[a, b], [c, d]].

    Rows are alleles (alt / ref), columns are strands (forward / reverse).
    """
    row1 = a + b
    col1 = a + c
    total = row1 + c + d
    if total == 0 or row1 == 0 or col1 == 0 or col1 == total:
        return 1.0

    def prob(x: int) -> float:
        return exp(
            _log_binom(col1, x) + _log_binom(total - col1, row1 - x) - _log_binom(total, row1)
        )

    p_obs = prob(a)
    lo = max(0, col1 - (total - row1))
    hi = min(col1, row1)
    cumulative = 0.0
    for x in range(lo, hi + 1):
        px = prob(x)
        if px <= p_obs * (1 + 1e-9):
            cumulative += px
    return min(1.0, cumulative)


def _adjacent_homopolymer(sequence: str, pos: int) -> bool:
    """True when a >=6 bp identical-base run overlaps pos or pos+1."""
    n = len(sequence)

    def run_length(center: int) -> int:
        base = sequence[center]
        left = center
        while (
            left - 1 >= 0
            and sequence[left - 1] == base
            and center - (left - 1) <= _HOMOPOLYMER_WINDOW
        ):
            left -= 1
        right = center
        while (
            right + 1 < n
            and sequence[right + 1] == base
            and (right + 1) - center <= _HOMOPOLYMER_WINDOW
        ):
            right += 1
        return right - left + 1

    for anchor in (pos, min(pos + 1, n - 1)):
        if (
            0 <= anchor < n
            and sequence[anchor] not in ("N", "")
            and run_length(anchor) >= _HOMOPOLYMER_RUN
        ):
            return True
    return False


# ---------------------------------------------------------------------------
# pileup walking
# ---------------------------------------------------------------------------


def call_heteroplasmy_from_bam(
    bam_path: Path,
    fasta_records: list[tuple[str, str]],
    min_mapq: int,
    min_baseq: int,
    min_depth: int,
    min_alt_count: int,
    min_alt_fraction: float,
) -> list[VariantCall]:
    import pysam  # local import keeps module importable without pysam

    references = {name: seq.upper() for name, seq in fasta_records}
    calls: list[VariantCall] = []

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for ref_name in bam.references:
            reference_seq = references.get(ref_name)
            if reference_seq is None:
                continue

            for column in bam.pileup(
                ref_name,
                stepper="samtools",
                min_baseq=min_baseq,
                min_mq=min_mapq,
            ):
                calls.extend(
                    _calls_from_counts(
                        ref_name,
                        column.reference_pos,
                        reference_seq,
                        *_collect_column_stats(column, reference_seq, min_baseq),
                        min_depth=min_depth,
                        min_alt_count=min_alt_count,
                        min_alt_fraction=min_alt_fraction,
                    )
                )

    calls.sort(key=lambda call: (call.contig, call.position, call.alt))
    return calls


def _collect_column_stats(column, reference_seq: str, min_baseq: int):
    """Single pass over a pileup column gathering every statistic we need."""
    pos = column.reference_pos
    ref_base = reference_seq[pos] if pos < len(reference_seq) else "N"
    n = len(reference_seq)

    depth = 0
    snp_counts: Counter[str] = Counter()
    strand_ref = [0, 0]  # forward, reverse
    snp_strand: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    snp_rel_positions: dict[str, list[float]] = defaultdict(list)
    indel_counts: Counter[tuple[str, int, str]] = Counter()
    indel_strand: dict[tuple[str, int, str], list[int]] = defaultdict(lambda: [0, 0])

    for pread in column.pileups:
        alignment = pread.alignment
        if (
            alignment.is_unmapped
            or alignment.is_secondary
            or alignment.is_supplementary
            or alignment.is_duplicate
        ):
            continue
        if alignment.mapping_quality < 1:  # hard floor; MAPQ filter applied in pileup()
            continue
        if pread.is_refskip:
            continue

        query_pos = pread.query_position
        if query_pos is None:
            continue

        qualities = alignment.query_qualities
        if qualities is not None and qualities[query_pos] < min_baseq:
            continue
        sequence = alignment.query_sequence or ""
        rev = alignment.is_reverse
        query_length = max(alignment.query_length or 1, 1)

        # ---- indels anchored at this position -----------------------------
        indel = pread.indel
        if indel > 0:
            inserted = sequence[query_pos + 1 : query_pos + 1 + indel]
            key = ("INS", len(inserted), inserted)
            indel_counts[key] += 1
            indel_strand[key][1 if rev else 0] += 1
            continue  # insertion reads do not vote for SNP/base stats here
        if indel < 0:
            deleted_len = -indel
            deleted = reference_seq[pos + 1 : pos + 1 + deleted_len]
            key = ("DEL", deleted_len, deleted)
            indel_counts[key] += 1
            indel_strand[key][1 if rev else 0] += 1
            continue

        if pread.is_del:
            continue  # placeholder D displayed by pileup stepper

        # ---- ordinary base -------------------------------------------------
        depth += 1
        base = sequence[query_pos : query_pos + 1] or "N"
        if base == ref_base:
            strand_ref[1 if rev else 0] += 1
            continue
        snp_counts[base] += 1
        snp_strand[base][1 if rev else 0] += 1
        snp_rel_positions[base].append(query_pos / query_length)

    # silence unused warnings for vars consumed indirectly
    _ = n
    return (
        ref_base,
        depth,
        snp_counts,
        tuple(strand_ref),
        snp_strand,
        snp_rel_positions,
        indel_counts,
        indel_strand,
    )


def _calls_from_counts(
    contig: str,
    position: int,
    reference_seq: str,
    ref_base: str,
    depth: int,
    snp_counts: Counter[str],
    strand_ref: tuple[int, int],
    snp_strand: dict[str, list[int]],
    snp_rel_positions: dict[str, list[float]],
    indel_counts: Counter[tuple[str, int, str]],
    indel_strand: dict[tuple[str, int, str], list[int]],
    *,
    min_depth: int,
    min_alt_count: int,
    min_alt_fraction: float,
) -> list[VariantCall]:
    """Build VariantCall objects from per-column statistics."""
    if depth < min_depth:
        return []

    def meets(count: int) -> bool:
        return count >= min_alt_count and (count / depth) >= min_alt_fraction

    calls: list[VariantCall] = []

    # ---- SNPs --------------------------------------------------------------
    for alt_base, count in snp_counts.most_common():
        if not meets(count):
            continue
        fwd, rev = snp_strand[alt_base]
        rel_positions = snp_rel_positions.get(alt_base, [])
        mean_rel_pos = sum(rel_positions) / len(rel_positions) if rel_positions else None
        calls.append(
            VariantCall(
                contig=contig,
                position=position,
                ref=ref_base,
                alt=alt_base,
                depth=depth,
                alt_count=count,
                alt_fraction=count / depth,
                variant_type="SNP",
                strand_ref_fwd=strand_ref[0],
                strand_ref_rev=strand_ref[1],
                strand_alt_fwd=fwd,
                strand_alt_rev=rev,
                strand_bias_pvalue=fisher_exact_2x2(fwd, rev, strand_ref[0], strand_ref[1]),
                homopolymer_context=_adjacent_homopolymer(reference_seq, position),
                mean_alt_read_position=round(mean_rel_pos, 4) if mean_rel_pos is not None else None,
            )
        )

    # ---- indels (VCF-style anchored alleles) --------------------------------
    anchor = reference_seq[position] if position < len(reference_seq) else "N"
    for (op, length, seq_fragment), count in indel_counts.most_common():
        if not meets(count):
            continue
        fwd, rev = indel_strand[(op, length, seq_fragment)]
        if op == "INS":
            ref_allele = anchor
            alt_allele = anchor + seq_fragment
        else:  # DEL
            ref_allele = anchor + seq_fragment
            alt_allele = anchor
        calls.append(
            VariantCall(
                contig=contig,
                position=position,
                ref=ref_allele,
                alt=alt_allele,
                depth=depth,
                alt_count=count,
                alt_fraction=count / depth,
                variant_type=op,
                strand_ref_fwd=strand_ref[0],
                strand_ref_rev=strand_ref[1],
                strand_alt_fwd=fwd,
                strand_alt_rev=rev,
                strand_bias_pvalue=fisher_exact_2x2(fwd, rev, strand_ref[0], strand_ref[1]),
                homopolymer_context=_adjacent_homopolymer(reference_seq, position),
                mean_alt_read_position=None,
            )
        )

    return calls
