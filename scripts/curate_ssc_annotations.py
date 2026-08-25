#!/usr/bin/env python3
"""Curate SSC coordinates for circular organelle assemblies.

Runs HELIOS inverted-repeat detection over every FASTA in an input directory
and writes a TSV of SSC/LSC/IR coordinates that can be fed back into the
pipeline via ``--ssc-start/--ssc-end`` when auto-detection is not possible at
runtime. Genomes without a credible IR pair (e.g., highly diverged or
reduced-IR plastomes like Durio) are listed with NA coordinates and must be
supplied manually.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from organelle_pipeline.parsers import read_fasta_records
from organelle_pipeline.repeats import detect_inverted_repeats


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", type=Path, help="Directory containing circular FASTA files")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output TSV (default: <input_dir>/ssc_annotations.tsv)",
    )
    parser.add_argument("--min-ir", type=int, default=10_000)
    parser.add_argument("--retry-min-ir", type=int, default=2_500)
    args = parser.parse_args()

    output = args.output or args.input_dir / "ssc_annotations.tsv"
    headers = [
        "sample",
        "contig",
        "contig_length",
        "ir_a_start",
        "ir_a_end",
        "ir_a_length",
        "ir_b_start",
        "ir_b_end",
        "ir_b_length",
        "ssc_start",
        "ssc_end",
        "ssc_length",
        "lsc_length",
        "source",
        "notes",
    ]

    rows: list[list[str]] = []
    fasta_paths = sorted(args.input_dir.glob("*.fasta")) + sorted(args.input_dir.glob("*.fa"))
    for path in fasta_paths:
        sample = path.stem.replace("_mq40_final_circular", "").replace("_final_circular", "")
        records = read_fasta_records(path)
        if len(records) != 1:
            rows.append(
                [
                    sample,
                    ";".join(name for name, _ in records),
                    str(len(records)),
                    *["NA"] * 11,
                    "multi_record",
                    "isomers taken directly from records by the pipeline",
                ]
            )
            continue
        name, sequence = records[0]
        n = len(sequence)
        detection = detect_inverted_repeats(sequence, min_ir_len=args.min_ir)
        if not detection.ok and any("only" in w for w in detection.warnings):
            detection = detect_inverted_repeats(sequence, min_ir_len=args.retry_min_ir)

        if detection.ok and detection.ir_a and detection.ir_b:
            ssc_len = (detection.ssc_end - detection.ssc_start) % n
            lsc_len = (detection.lsc_end - detection.lsc_start) % n
            total = detection.ir_a.length + detection.ir_b.length + ssc_len + lsc_len
            notes = "coordinates sum to contig length" if total == n else f"SUM MISMATCH ({total})"
            if detection.warnings:
                notes += "; " + "; ".join(detection.warnings)
            rows.append(
                [
                    sample,
                    name,
                    str(n),
                    str(detection.ir_a.start),
                    str(detection.ir_a.end),
                    str(detection.ir_a.length),
                    str(detection.ir_b.start),
                    str(detection.ir_b.end),
                    str(detection.ir_b.length),
                    str(detection.ssc_start),
                    str(detection.ssc_end),
                    str(ssc_len),
                    str(lsc_len),
                    "auto_detected",
                    notes,
                ]
            )
        else:
            rows.append(
                [
                    sample,
                    name,
                    str(n),
                    *["NA"] * 10,
                    "undetermined",
                    "no credible IR pair detected; supply --ssc-start/--ssc-end manually; "
                    + "; ".join(detection.warnings),
                ]
            )

    output.write_text("\t".join(headers) + "\n" + "".join("\t".join(r) + "\n" for r in rows))
    print(f"wrote {len(rows)} entries -> {output}")
    for row in rows:
        print(f"  {row[0]:12s} SSC={row[9]}-{row[10]} ({row[11]}) [{row[13]}]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
