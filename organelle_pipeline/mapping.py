import shutil
from pathlib import Path

from organelle_pipeline.utils import run_command


def choose_aligner(preferred: str, average_read_length: float) -> str:
    if preferred in {"bwa", "minimap2"}:
        return preferred
    if preferred != "auto":
        raise ValueError("Aligner must be one of: auto, bwa, minimap2")
    return "minimap2" if average_read_length >= 800 else "bwa"


def run_mapping(
    aligner: str,
    fasta_path: Path,
    fastq_files: list[Path],
    output_bam: Path,
    threads: int,
) -> Path:
    output_bam.parent.mkdir(parents=True, exist_ok=True)

    if aligner == "bwa":
        _require_tools(["bwa", "samtools"])
        _build_bwa_index(fasta_path)
        temp_sam = output_bam.with_suffix(".tmp.sam")
        command = ["bwa", "mem", "-t", str(threads), str(fasta_path)] + [str(p) for p in fastq_files]
        _run_to_file(command, temp_sam)
        run_command(["samtools", "sort", "-@", str(threads), "-o", str(output_bam), str(temp_sam)])
        temp_sam.unlink(missing_ok=True)
    elif aligner == "minimap2":
        _require_tools(["minimap2", "samtools"])
        index_path = fasta_path.with_suffix(".mmi")
        if not index_path.exists():
            run_command(["minimap2", "-d", str(index_path), str(fasta_path)])
        preset = "map-ont" if _looks_long_read(fastq_files) else "sr"
        temp_sam = output_bam.with_suffix(".tmp.sam")
        command = [
            "minimap2",
            "-a",
            "-x",
            preset,
            "-t",
            str(threads),
            str(index_path),
        ] + [str(p) for p in fastq_files]
        _run_to_file(command, temp_sam)
        run_command(["samtools", "sort", "-@", str(threads), "-o", str(output_bam), str(temp_sam)])
        temp_sam.unlink(missing_ok=True)
    else:
        raise ValueError(f"Unsupported aligner: {aligner}")

    run_command(["samtools", "index", str(output_bam)])
    return output_bam


def _looks_long_read(fastq_files: list[Path], probe_reads: int = 100) -> bool:
    from organelle_pipeline.parsers import read_fastq_sequences

    total = 0
    count = 0
    for sequence in read_fastq_sequences(fastq_files, limit=probe_reads):
        total += len(sequence)
        count += 1
    if count == 0:
        return False
    return (total / count) >= 800


def _require_tools(tools: list[str]) -> None:
    missing = [tool for tool in tools if shutil.which(tool) is None]
    if missing:
        joined = ", ".join(missing)
        raise RuntimeError(f"Missing required external tools: {joined}")


def _build_bwa_index(fasta_path: Path) -> None:
    bwt_path = fasta_path.with_suffix(fasta_path.suffix + ".bwt")
    if bwt_path.exists():
        return
    run_command(["bwa", "index", str(fasta_path)])


def _run_to_file(command: list[str], output_path: Path) -> None:
    import subprocess

    with output_path.open("w", encoding="utf-8") as handle:
        completed = subprocess.run(command, stdout=handle, stderr=subprocess.PIPE, text=True)
    if completed.returncode != 0:
        stderr = completed.stderr.strip() if completed.stderr else "unknown error"
        raise RuntimeError(f"Command failed: {' '.join(command)}\n{stderr}")
