from collections import Counter
from pathlib import Path

from organelle_pipeline.models import VariantCall


def call_heteroplasmy_from_bam(
    bam_path: Path,
    fasta_records: list[tuple[str, str]],
    min_mapq: int,
    min_baseq: int,
    min_depth: int,
    min_alt_count: int,
    min_alt_fraction: float,
) -> list[VariantCall]:
    try:
        pysam = __import__("pysam")
    except ImportError as exc:
        raise RuntimeError("pysam is required for heteroplasmy calling") from exc

    reference_map = {name: sequence for name, sequence in fasta_records}
    calls: list[VariantCall] = []

    with pysam.AlignmentFile(str(bam_path), "rb") as bam_file:
        for contig in bam_file.references:
            if contig not in reference_map:
                continue
            reference_sequence = reference_map[contig]
            pileup_iter = bam_file.pileup(
                contig,
                truncate=True,
                min_base_quality=min_baseq,
                min_mapping_quality=min_mapq,
                stepper="samtools",
            )
            for column in pileup_iter:
                position_1 = column.reference_pos + 1
                if position_1 > len(reference_sequence):
                    continue
                ref_base = reference_sequence[position_1 - 1]
                base_counter: Counter[str] = Counter()
                indel_counter: Counter[str] = Counter()
                depth = 0

                for pileup_read in column.pileups:
                    alignment = pileup_read.alignment
                    if alignment.mapping_quality < min_mapq:
                        continue

                    if pileup_read.is_del:
                        indel_counter["DEL"] += 1
                        depth += 1
                        continue
                    if pileup_read.is_refskip:
                        continue

                    query_pos = pileup_read.query_position
                    if query_pos is None:
                        continue

                    base = alignment.query_sequence[query_pos].upper()
                    quality = alignment.query_qualities[query_pos]
                    if quality < min_baseq:
                        continue

                    depth += 1
                    base_counter[base] += 1
                    if pileup_read.indel > 0:
                        indel_counter["INS"] += 1

                calls.extend(
                    _calls_from_counts(
                        contig=contig,
                        position=position_1,
                        ref_base=ref_base,
                        base_counter=base_counter,
                        indel_counter=indel_counter,
                        depth=depth,
                        min_depth=min_depth,
                        min_alt_count=min_alt_count,
                        min_alt_fraction=min_alt_fraction,
                    )
                )

    calls.sort(key=lambda item: (item.contig, item.position, item.alt))
    return calls


def _calls_from_counts(
    contig: str,
    position: int,
    ref_base: str,
    base_counter: Counter[str],
    indel_counter: Counter[str],
    depth: int,
    min_depth: int,
    min_alt_count: int,
    min_alt_fraction: float,
) -> list[VariantCall]:
    if depth < min_depth:
        return []

    out: list[VariantCall] = []
    for base, count in base_counter.items():
        if base == ref_base:
            continue
        fraction = count / depth
        if count < min_alt_count or fraction < min_alt_fraction:
            continue
        out.append(
            VariantCall(
                contig=contig,
                position=position,
                ref=ref_base,
                alt=base,
                depth=depth,
                alt_count=count,
                alt_fraction=fraction,
                variant_type="SNP",
            )
        )

    for indel_type in ("INS", "DEL"):
        count = indel_counter.get(indel_type, 0)
        fraction = count / depth if depth else 0.0
        if count < min_alt_count or fraction < min_alt_fraction:
            continue
        out.append(
            VariantCall(
                contig=contig,
                position=position,
                ref=ref_base,
                alt=indel_type,
                depth=depth,
                alt_count=count,
                alt_fraction=fraction,
                variant_type="INDEL",
            )
        )
    return out
