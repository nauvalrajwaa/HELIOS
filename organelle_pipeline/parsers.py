import re
from pathlib import Path

from organelle_pipeline.utils import open_text


def read_fasta_records(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name = ""
    chunks: list[str] = []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    records.append((name, "".join(chunks).upper()))
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name:
        records.append((name, "".join(chunks).upper()))
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records


def read_fastq_sequences(paths: list[Path], limit: int | None = None):
    emitted = 0
    for path in paths:
        with open_text(path) as handle:
            while True:
                header = handle.readline()
                if not header:
                    break
                sequence = handle.readline().strip()
                plus = handle.readline()
                quality = handle.readline()
                if not plus or not quality:
                    raise ValueError(f"Malformed FASTQ record in {path}")
                yield sequence.upper()
                emitted += 1
                if limit is not None and emitted >= limit:
                    return


def estimate_read_length(paths: list[Path], limit: int = 2000) -> float:
    total = 0
    count = 0
    for sequence in read_fastq_sequences(paths, limit=limit):
        total += len(sequence)
        count += 1
    if count == 0:
        raise ValueError("No FASTQ reads found")
    return total / count


def infer_ssc_region(annotation_path: Path | None) -> tuple[int, int] | None:
    if annotation_path is None:
        return None
    suffix = annotation_path.suffix.lower()
    if suffix in {".gff", ".gff3"}:
        return _infer_from_gff(annotation_path)
    if suffix in {".gb", ".gbk", ".genbank"}:
        return _infer_from_gbk(annotation_path)
    return None


def _infer_from_gff(path: Path) -> tuple[int, int] | None:
    with open_text(path) as handle:
        for raw in handle:
            if not raw or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature_type = fields[2].lower()
            attrs = fields[8].lower()
            keyspace = f"{feature_type} {attrs}"
            if "ssc" in keyspace or "small_single_copy" in keyspace:
                start = int(fields[3])
                end = int(fields[4])
                if start > end:
                    start, end = end, start
                return start, end
    return None


def _infer_from_gbk(path: Path) -> tuple[int, int] | None:
    pattern = re.compile(r"(\d+)\.\.(\d+)")
    with open_text(path) as handle:
        for raw in handle:
            lower = raw.lower()
            if "ssc" not in lower and "small_single_copy" not in lower:
                continue
            match = pattern.search(raw)
            if match:
                start = int(match.group(1))
                end = int(match.group(2))
                if start > end:
                    start, end = end, start
                return start, end
    return None
