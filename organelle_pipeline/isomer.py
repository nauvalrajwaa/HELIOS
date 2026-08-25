"""Isomer construction and read-backed quantification.

Two complementary estimators are provided:

``quantify_isomers_by_remapping``
    Primary estimator. Reads are mapped against BOTH isomer references; the
    orientation of alignments falling inside the SSC core discriminates which
    structural isoform each read originated from (a read whose molecule
    carries the flipped SSC aligns antisense to the A reference inside the
    SSC, and sense to the B reference). Junction-spanning reads are excluded
    via a buffer so only unambiguous molecules vote.

``quantify_isomers_from_fastq``
    Secondary cross-check. Unique k-mer voting between the two isomer
    sequences (kept from HELIOS v0.1).
"""

from __future__ import annotations

from dataclasses import dataclass, field

from organelle_pipeline.models import IsomerQuantResult
from organelle_pipeline.parsers import read_fastq_sequences
from organelle_pipeline.repeats import IrDetection, invert_segment
from organelle_pipeline.utils import reverse_complement


@dataclass(slots=True)
class IsomerPair:
    """Isomer A/B sequences plus provenance of the SSC definition."""

    name_a: str
    sequence_a: str
    name_b: str
    sequence_b: str
    ssc_start: int | None = None  # 0-based inclusive on the A frame (may wrap)
    ssc_end: int | None = None  # 0-based exclusive on the A frame (may wrap)
    ssc_source: str = "unknown"  # user_override | annotation | auto_detected | two_record
    assumptions: list[str] = field(default_factory=list)
    contig_name: str | None = None  # original record name for dual-reference remapping


def build_isomer_candidates(
    fasta_records: list[tuple[str, str]],
    ssc_region: tuple[int, int] | None = None,
    ssc_source: str = "annotation",
    detection: IrDetection | None = None,
) -> IsomerPair:
    """Construct the isomer pair from an assembly.

    Priority: explicit ``ssc_region`` (CLI override / curated annotation) >
    auto-detected IR boundaries (``detection``) > two-record FASTA shortcut.
    Raises ``ValueError`` with actionable guidance when nothing applies.
    """
    if not fasta_records:
        raise ValueError("no FASTA records supplied")

    if len(fasta_records) >= 2:
        (name_a, seq_a), (name_b, seq_b) = fasta_records[0], fasta_records[1]
        return IsomerPair(
            name_a=name_a,
            sequence_a=seq_a,
            name_b=name_b,
            sequence_b=seq_b,
            ssc_start=None,
            ssc_end=None,
            ssc_source="two_record",
            assumptions=[
                "Two or more FASTA records detected; first two treated as Isomer A and B.",
                "Orientation-based remapping unavailable without SSC coordinates; "
                "k-mer voting used instead.",
            ],
        )

    name, sequence = fasta_records[0]
    n = len(sequence)

    # ---- explicit region (user override or curated annotation) -------------
    if ssc_region:
        start, end = ssc_region
        start = max(0, int(start))
        end = min(n, int(end))
        if end - start < max(200, n // 40):
            raise ValueError(
                f"SSC region {start}-{end} too small ({end - start} bp) for contig of {n} bp"
            )
        isomer_b_sequence = invert_segment(sequence, start, end)
        return IsomerPair(
            name_a=f"{name}_A",
            contig_name=name,
            sequence_a=sequence,
            name_b=f"{name}_B",
            sequence_b=isomer_b_sequence,
            ssc_start=start,
            ssc_end=end,
            ssc_source=ssc_source,
            assumptions=[
                f"SSC region {start}-{end} ({ssc_source}) inverted to build Isomer B.",
            ],
        )

    # ---- auto-detected inverted repeat boundaries --------------------------
    if detection is not None and detection.ok and detection.ssc_start is not None:
        start, end = detection.ssc_start, detection.ssc_end
        isomer_b_sequence = invert_segment(sequence, start, end)
        notes = [
            f"SSC region [{start}, {end}) derived from self-detected inverted repeats "
            "(IRA/IRB); Isomer B built by inverting this arc."
        ]
        for warning in detection.warnings:
            notes.append(f"IR detection warning: {warning}")
        return IsomerPair(
            name_a=f"{name}_A",
            contig_name=name,
            sequence_a=sequence,
            name_b=f"{name}_B",
            sequence_b=isomer_b_sequence,
            ssc_start=start,
            ssc_end=end,
            ssc_source="auto_detected",
            assumptions=notes,
        )

    raise ValueError(
        "Could not determine SSC boundaries for a single-record assembly "
        f"({name}, {n} bp). No inverted-repeat pair was detectable. Provide "
        "--ssc-start/--ssc-end explicitly or a GFF3/GenBank annotation carrying the "
        "SSC feature."
    )


# ---------------------------------------------------------------------------
# primary estimator: dual-reference remapping
# ---------------------------------------------------------------------------


@dataclass(slots=True)
class RemapQuantResult:
    fraction_a: float
    fraction_b: float
    sense_a: int  # reads supporting Isomer A when aligned to reference A
    antisense_a: int  # reads supporting Isomer B when aligned to reference A
    total_ssc_reads_a: int
    sense_b: int  # reads supporting Isomer B when aligned to reference B
    antisense_b: int  # reads supporting Isomer A when aligned to reference B
    total_ssc_reads_b: int
    agreement_delta: float  # |f_B(refA) - f_B(refB)|; small = high confidence
    junction_buffer: int
    min_mapq: int


def _ssc_core_pieces(
    ssc_start: int, ssc_end: int, n: int, junction_buffer: int
) -> list[tuple[int, int]]:
    """Linear pieces of the SSC interior after removing junction buffers."""
    core_len = (ssc_end - ssc_start) % n
    lo = ssc_start + junction_buffer
    hi = ssc_start + core_len - junction_buffer  # may exceed n (fine)
    if hi - lo <= 0:
        return []
    pieces: list[tuple[int, int]] = []
    pos = lo
    while pos < hi:
        piece_end = min(((pos // n) + 1) * n, hi)
        pieces.append((pos % n, piece_end % n or n))
        pos = piece_end
    return pieces


def _count_orientations(bam_path, ref_name: str, regions, min_mapq: int) -> tuple[int, int]:
    import pysam  # local import keeps module importable without pysam

    sense = antisense = 0
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for region_start, region_end in regions:
            for read in bam.fetch(ref_name, region_start, region_end):
                if (
                    read.is_unmapped
                    or read.is_secondary
                    or read.is_supplementary
                    or read.is_duplicate
                ):
                    continue
                if read.mapping_quality < min_mapq:
                    continue
                # keep only reads fully contained in the SSC core
                if read.reference_start < region_start or read.reference_end > region_end:
                    continue
                if read.is_paired and not read.is_read1:
                    continue  # count MOLECULES via first-in-pair only;
                # every proper FR pair otherwise contributes exactly one
                # forward and one reverse mate and cancels itself out.
                if read.is_reverse:
                    antisense += 1
                else:
                    sense += 1
    return sense, antisense


def quantify_isomers_by_remapping(
    bam_a_path,
    bam_b_path,
    ref_name: str,
    ref_length: int,
    ssc_start: int,
    ssc_end: int,
    *,
    junction_buffer: int = 250,
    min_mapq: int = 20,
) -> RemapQuantResult | None:
    """Estimate isomer proportions from alignment orientation inside the SSC.

    Against reference A: sense reads = Isomer A molecules, antisense = Isomer B.
    Against reference B (same coordinate frame): sense = B, antisense = A.
    Both frames must agree; their disagreement is reported as QC.
    """
    regions = _ssc_core_pieces(ssc_start, ssc_end, ref_length, junction_buffer)
    if not regions:
        return None

    sense_a, antisense_a = _count_orientations(bam_a_path, ref_name, regions, min_mapq)
    sense_b, antisense_b = _count_orientations(bam_b_path, ref_name, regions, min_mapq)

    informative_a = sense_a + antisense_a
    informative_b = sense_b + antisense_b
    if informative_a == 0 or informative_b == 0:
        return None

    frac_b_from_a = antisense_a / informative_a
    frac_b_from_b = sense_b / informative_b
    frac_b = (frac_b_from_a + frac_b_from_b) / 2.0
    return RemapQuantResult(
        fraction_a=1.0 - frac_b,
        fraction_b=frac_b,
        sense_a=sense_a,
        antisense_a=antisense_a,
        total_ssc_reads_a=informative_a,
        sense_b=sense_b,
        antisense_b=antisense_b,
        total_ssc_reads_b=informative_b,
        agreement_delta=abs(frac_b_from_a - frac_b_from_b),
        junction_buffer=junction_buffer,
        min_mapq=min_mapq,
    )


# ---------------------------------------------------------------------------
# secondary estimator: unique k-mer voting
# ---------------------------------------------------------------------------


def _build_unique_kmers(seq_a: str, seq_b: str, kmer_size: int) -> tuple[set[str], set[str]]:
    kmers_a = {seq_a[i : i + kmer_size] for i in range(len(seq_a) - kmer_size + 1)}
    kmers_b = {seq_b[i : i + kmer_size] for i in range(len(seq_b) - kmer_size + 1)}
    return kmers_a - kmers_b, kmers_b - kmers_a


def _count_hits(read: str, kmers: set[str], kmer_size: int) -> int:
    return sum(1 for i in range(len(read) - kmer_size + 1) if read[i : i + kmer_size] in kmers)


def quantify_isomers_from_fastq(
    fastq_files,
    isomer_pair: IsomerPair,
    kmer_size: int = 31,
    min_hits: int = 2,
    read_limit: int = 250_000,
    assumptions: list[str] | None = None,
) -> IsomerQuantResult:
    """Unique k-mer voting estimator (secondary cross-check)."""
    unique_a, unique_b = _build_unique_kmers(
        isomer_pair.sequence_a.upper(), isomer_pair.sequence_b.upper(), kmer_size
    )
    assigned_a = assigned_b = ambiguous = unassigned = seen = 0
    for record in read_fastq_sequences(fastq_files, limit=read_limit):
        seen += 1
        read = record.upper()
        hits_a_fwd = _count_hits(read, unique_a, kmer_size)
        hits_b_fwd = _count_hits(read, unique_b, kmer_size)
        rc_read = reverse_complement(read)
        hits_a_rc = _count_hits(rc_read, unique_a, kmer_size)
        hits_b_rc = _count_hits(rc_read, unique_b, kmer_size)
        # orientation-first: read the sequence in whichever direction carries
        # more total support, THEN compare isomers within that single frame.
        if hits_a_rc + hits_b_rc > hits_a_fwd + hits_b_fwd:
            hits_a, hits_b = hits_a_rc, hits_b_rc
        else:
            hits_a, hits_b = hits_a_fwd, hits_b_fwd
        if max(hits_a, hits_b) < min_hits:
            unassigned += 1
        elif hits_a > hits_b:
            assigned_a += 1
        elif hits_b > hits_a:
            assigned_b += 1
        else:
            ambiguous += 1

    informative = assigned_a + assigned_b
    fraction_a = assigned_a / informative if informative else 0.0
    fraction_b = assigned_b / informative if informative else 0.0

    notes = list(assumptions or [])
    notes.append(
        "K-mer voting uses unique k-mers (symmetric difference) at "
        f"k={kmer_size}; orientation-aware."
    )
    return IsomerQuantResult(
        isomer_a_name=isomer_pair.name_a,
        isomer_b_name=isomer_pair.name_b,
        assigned_a=assigned_a,
        assigned_b=assigned_b,
        ambiguous=ambiguous,
        unassigned=unassigned,
        total_reads_seen=seen,
        isomer_a_fraction=round(fraction_a, 6),
        isomer_b_fraction=round(fraction_b, 6),
        method="unique_kmer_voting",
        assumptions=notes,
    )


def write_isomer_gfa(
    path,
    isomer_pair: IsomerPair,
    quant_result: IsomerQuantResult | None,
) -> None:
    """Write a GFA2-style graph with REAL segment sequences.

    Segments carry the actual isomer sequences (renderable by Bandage/odgi);
    RC:i holds assigned-read counts and PR:f the proportion custom tags.
    """
    path.parent.mkdir(parents=True, exist_ok=True)

    def _segment_line(seg_name: str, sequence: str, reads: int, proportion: float) -> str:
        return f"S\t{seg_name}\t{sequence}\tRC:i:{reads}\tPR:f:{proportion:.6f}\n"

    lines = ["H\tVN:Z:2.0\n"]
    if quant_result is not None:
        lines.append(
            _segment_line(
                isomer_pair.name_a,
                isomer_pair.sequence_a,
                quant_result.assigned_a,
                quant_result.isomer_a_fraction,
            )
        )
        lines.append(
            _segment_line(
                isomer_pair.name_b,
                isomer_pair.sequence_b,
                quant_result.assigned_b,
                quant_result.isomer_b_fraction,
            )
        )
        lines.append(f"L\t{isomer_pair.name_a}\t+\t{isomer_pair.name_b}\t-\t*\n")
    else:
        lines.append(_segment_line(isomer_pair.name_a, isomer_pair.sequence_a, 0, float("nan")))
        lines.append(_segment_line(isomer_pair.name_b, isomer_pair.sequence_b, 0, float("nan")))
    path.write_text("".join(lines))
