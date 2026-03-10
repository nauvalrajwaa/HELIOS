from pathlib import Path

from organelle_pipeline.models import IsomerQuantResult
from organelle_pipeline.parsers import read_fastq_sequences
from organelle_pipeline.utils import reverse_complement


def build_isomer_candidates(
    fasta_records: list[tuple[str, str]],
    ssc_region: tuple[int, int] | None,
) -> tuple[tuple[str, str], tuple[str, str], list[str]]:
    if len(fasta_records) >= 2:
        assumptions = [
            "Two or more FASTA records detected; first two records treated as Isomer A and Isomer B candidates.",
            "Isomer quantification uses record-specific unique k-mers.",
        ]
        return fasta_records[0], fasta_records[1], assumptions

    name, sequence = fasta_records[0]
    if ssc_region:
        start_1, end_1 = ssc_region
        start = max(0, start_1 - 1)
        end = min(len(sequence), end_1)
        if start < end:
            isomer_b_sequence = _invert_segment(sequence, start, end)
            assumptions = [
                f"Single FASTA record detected; inferred SSC region {start_1}-{end_1} from annotation and created Isomer B by inversion.",
                "Isomer quantification uses unique k-mers from Isomer A and inferred Isomer B.",
            ]
            return (f"{name}_A", sequence), (f"{name}_B", isomer_b_sequence), assumptions

    segment_start = int(len(sequence) * 0.35)
    segment_end = int(len(sequence) * 0.65)
    isomer_b_sequence = _invert_segment(sequence, segment_start, segment_end)
    assumptions = [
        "Single FASTA record detected and SSC annotation not found; Isomer B approximated by inverting the middle 30% segment.",
        "This fallback supports comparative structural quantification but should be replaced by curated isomer references when available.",
    ]
    return (f"{name}_A", sequence), (f"{name}_B", isomer_b_sequence), assumptions


def quantify_isomers_from_fastq(
    fastq_files: list[Path],
    isomer_a: tuple[str, str],
    isomer_b: tuple[str, str],
    kmer_size: int,
    min_hits: int,
    read_limit: int,
    assumptions: list[str],
) -> IsomerQuantResult:
    isomer_a_name, seq_a = isomer_a
    isomer_b_name, seq_b = isomer_b

    unique_a, unique_b = _build_unique_kmers(seq_a, seq_b, kmer_size)

    assigned_a = 0
    assigned_b = 0
    ambiguous = 0
    unassigned = 0
    total = 0

    for read_sequence in read_fastq_sequences(fastq_files, limit=read_limit):
        total += 1
        a_forward = _count_hits(read_sequence, unique_a, kmer_size)
        b_forward = _count_hits(read_sequence, unique_b, kmer_size)
        read_rc = reverse_complement(read_sequence)
        a_reverse = _count_hits(read_rc, unique_a, kmer_size)
        b_reverse = _count_hits(read_rc, unique_b, kmer_size)

        a_hits, b_hits = _select_orientation((a_forward, b_forward), (a_reverse, b_reverse))

        if a_hits < min_hits and b_hits < min_hits:
            unassigned += 1
            continue
        if a_hits == b_hits:
            ambiguous += 1
            continue
        if a_hits > b_hits:
            assigned_a += 1
        else:
            assigned_b += 1

    informative = assigned_a + assigned_b
    if informative == 0:
        frac_a = 0.0
        frac_b = 0.0
    else:
        frac_a = assigned_a / informative
        frac_b = assigned_b / informative

    return IsomerQuantResult(
        isomer_a_name=isomer_a_name,
        isomer_b_name=isomer_b_name,
        assigned_a=assigned_a,
        assigned_b=assigned_b,
        ambiguous=ambiguous,
        unassigned=unassigned,
        total_reads_seen=total,
        isomer_a_fraction=frac_a,
        isomer_b_fraction=frac_b,
        method="unique_kmer_voting",
        assumptions=assumptions,
    )


def write_isomer_gfa(path: Path, isomer_result: IsomerQuantResult) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    a_pct = round(isomer_result.isomer_a_fraction * 100.0, 4)
    b_pct = round(isomer_result.isomer_b_fraction * 100.0, 4)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("H\tVN:Z:1.0\n")
        handle.write(
            f"S\t{isomer_result.isomer_a_name}\t*\tRC:i:{isomer_result.assigned_a}\tPR:f:{a_pct}\n"
        )
        handle.write(
            f"S\t{isomer_result.isomer_b_name}\t*\tRC:i:{isomer_result.assigned_b}\tPR:f:{b_pct}\n"
        )
        handle.write(
            f"L\t{isomer_result.isomer_a_name}\t+\t{isomer_result.isomer_b_name}\t+\t0M\n"
        )
        handle.write(
            f"L\t{isomer_result.isomer_b_name}\t+\t{isomer_result.isomer_a_name}\t+\t0M\n"
        )
        handle.write(
            f"P\tIsomerPath\t{isomer_result.isomer_a_name}+,{isomer_result.isomer_b_name}+\t*\n"
        )


def _invert_segment(sequence: str, start: int, end: int) -> str:
    return sequence[:start] + reverse_complement(sequence[start:end]) + sequence[end:]


def _build_unique_kmers(seq_a: str, seq_b: str, kmer_size: int) -> tuple[set[str], set[str]]:
    kmers_a = _kmers(seq_a, kmer_size)
    kmers_b = _kmers(seq_b, kmer_size)
    unique_a = kmers_a - kmers_b
    unique_b = kmers_b - kmers_a
    return unique_a, unique_b


def _kmers(sequence: str, kmer_size: int) -> set[str]:
    if kmer_size <= 0:
        raise ValueError("kmer_size must be > 0")
    if len(sequence) < kmer_size:
        return set()
    out: set[str] = set()
    for index in range(0, len(sequence) - kmer_size + 1):
        kmer = sequence[index : index + kmer_size]
        if "N" in kmer:
            continue
        out.add(kmer)
    return out


def _count_hits(sequence: str, kmer_set: set[str], kmer_size: int) -> int:
    if len(sequence) < kmer_size or not kmer_set:
        return 0
    hits = 0
    for index in range(0, len(sequence) - kmer_size + 1):
        if sequence[index : index + kmer_size] in kmer_set:
            hits += 1
    return hits


def _select_orientation(forward: tuple[int, int], reverse: tuple[int, int]) -> tuple[int, int]:
    f_a, f_b = forward
    r_a, r_b = reverse
    f_margin = abs(f_a - f_b)
    r_margin = abs(r_a - r_b)
    if r_margin > f_margin:
        return reverse
    if f_margin > r_margin:
        return forward
    if (r_a + r_b) > (f_a + f_b):
        return reverse
    return forward
