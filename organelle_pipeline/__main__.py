import argparse
from pathlib import Path

from organelle_pipeline.models import PipelineConfig
from organelle_pipeline.pipeline import run_pipeline


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="helios",
        description="HELIOS: HEteroplasmy and Isomer Locator for Organelle Sequences",
    )
    parser.add_argument("--fasta", type=Path, required=True, help="Circular organelle FASTA")
    parser.add_argument("--annotation", type=Path, default=None, help="Optional GFF3 or GBK annotation")
    parser.add_argument(
        "--fastq",
        type=Path,
        nargs="+",
        default=None,
        help="Optional raw FASTQ file(s) for read-backed isomer quantification and heteroplasmy calling",
    )
    parser.add_argument("--output", type=Path, required=True, help="Output directory")
    parser.add_argument("--sample-name", type=str, default="sample", help="Sample name for output naming")
    parser.add_argument("--aligner", choices=["auto", "bwa", "minimap2"], default="auto")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--min-baseq", type=int, default=20)
    parser.add_argument("--min-depth", type=int, default=200)
    parser.add_argument("--min-alt-count", type=int, default=10)
    parser.add_argument("--min-alt-fraction", type=float, default=0.01)
    parser.add_argument("--kmer-size", type=int, default=31)
    parser.add_argument("--min-isomer-hits", type=int, default=2)
    parser.add_argument("--skip-mapping", action="store_true", help="Skip BAM remapping and heteroplasmy calling")
    parser.add_argument(
        "--read-limit-for-stats",
        type=int,
        default=2000,
        help="Maximum reads used for read-length/aligner inference",
    )
    parser.add_argument(
        "--read-limit-for-isomer",
        type=int,
        default=250000,
        help="Maximum reads processed for isomer quantification",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    config = PipelineConfig(
        fasta=args.fasta,
        annotation=args.annotation,
        fastq_files=args.fastq or [],
        output_dir=args.output,
        sample_name=args.sample_name,
        aligner=args.aligner,
        threads=args.threads,
        min_mapq=args.min_mapq,
        min_baseq=args.min_baseq,
        min_depth=args.min_depth,
        min_alt_count=args.min_alt_count,
        min_alt_fraction=args.min_alt_fraction,
        kmer_size=args.kmer_size,
        min_isomer_hits=args.min_isomer_hits,
        skip_mapping=args.skip_mapping,
        read_limit_for_stats=args.read_limit_for_stats,
        read_limit_for_isomer=args.read_limit_for_isomer,
    )

    result = run_pipeline(config)
    print(f"Completed sample: {result.sample_name}")
    print(f"Report: {result.html_report}")
    print(f"GFA: {result.gfa_graph}")
    print(f"Summary: {result.summary_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
