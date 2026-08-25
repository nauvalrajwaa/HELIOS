"""End-to-end pipeline tests.

The integration test builds a synthetic plastome with engineered inverted
repeats, simulates paired-end reads containing a KNOWN mixture of Isomer A
and Isomer B molecules inside the SSC, runs the full HELIOS pipeline
(auto-detection -> dual-reference remapping -> heteroplasmy), and asserts the
ground-truth proportions are recovered.
"""

import json
import random
import shutil

import pytest

from organelle_pipeline.models import PipelineConfig
from organelle_pipeline.pipeline import run_pipeline
from organelle_pipeline.utils import reverse_complement

requires_aligners = pytest.mark.skipif(
    shutil.which("bwa") is None or shutil.which("samtools") is None,
    reason="bwa/samtools must be installed for end-to-end tests",
)

LSC_LEN = 20_000
IR_LEN = 8_000
SSC_LEN = 6_000
READ_LEN = 150
INSERT = 300


def _build_synthetic_genome():
    random.seed(5)
    ir = "".join(random.choice("ACGT") for _ in range(IR_LEN))
    ssc = "".join(random.choice("ACGT") for _ in range(SSC_LEN))
    lsc = "".join(random.choice("ACGT") for _ in range(LSC_LEN))
    return lsc + ir + ssc + reverse_complement(ir)


def _write_fasta(path, name, sequence):
    path.write_text(f">{name}\n{sequence}\n")


def _pair_from_molecule(molecule, start):
    end = start + INSERT
    if end > len(molecule):
        return None
    r1 = molecule[start : start + READ_LEN]
    r2 = reverse_complement(molecule[end - READ_LEN : end])
    if "N" in r1 or "N" in r2:
        return None
    return r1, r2


def _simulate_reads(
    genome,
    ssc_start,
    ssc_end,
    n_a_ssc=700,
    n_b_ssc=300,
    n_background=1500,
    seed=17,
):
    """Return list of (r1, r2) pairs: known 30%% B-molecules inside the SSC."""
    rng = random.Random(seed)
    pairs = []

    # informative A molecules: sense orientation against reference A
    lo = ssc_start + 200
    hi = ssc_end - INSERT - 200
    while len(pairs) < n_a_ssc:
        pair = _pair_from_molecule(genome, rng.randint(lo, hi))
        if pair:
            pairs.append(pair)

    # informative B molecules: their SSC equals the reverse complement arc
    b_molecule = reverse_complement(genome[ssc_start:ssc_end])
    made_b = 0
    while made_b < n_b_ssc:
        pair = _pair_from_molecule(b_molecule, rng.randint(lo - ssc_start, hi - ssc_start))
        if pair:
            pairs.append(pair)
            made_b += 1

    # background reads outside the SSC (uninformative for orientation)
    made_bg = 0
    while made_bg < n_background:
        if rng.random() < 0.5:
            pos = rng.randint(0, ssc_start - INSERT - 1)
        else:
            pos = rng.randint(ssc_end, len(genome) - INSERT - 1)
        pair = _pair_from_molecule(genome, pos)
        if pair:
            pairs.append(pair)
            made_bg += 1

    rng.shuffle(pairs)
    return pairs


def _write_fastq(path, pairs, first=True):
    suffix = "/1" if first else "/2"
    lines = []
    idx_name = 0 if first else 1
    for i, pair in enumerate(pairs):
        lines += [
            f"@synthetic_{i}{suffix}",
            pair[idx_name],
            "+",
            "I" * len(pair[idx_name]),
        ]
    path.write_text("\n".join(lines) + "\n")


def _make_config(tmp_path, fasta, fastqs, **overrides):
    defaults = {
        "fasta": fasta,
        "annotation": None,
        "fastq_files": fastqs,
        "output_dir": tmp_path / "out",
        "sample_name": "synthetic",
        "aligner": "bwa",
        "threads": 2,
        "min_mapq": 20,
        "min_baseq": 20,
        "min_depth": 10,
        "min_alt_count": 3,
        "min_alt_fraction": 0.01,
        "kmer_size": 21,
        "min_isomer_hits": 2,
        "method": "both",
        "ssc_start": None,
        "ssc_end": None,
        "junction_buffer": 100,
        "read_limit_for_stats": 500,
        "read_limit_for_isomer": 100_000,
    }
    defaults.update(overrides)
    return PipelineConfig(**defaults)


def test_missing_fastq_raises_actionable_error(tmp_path):
    fasta = tmp_path / "genome.fasta"
    _write_fasta(fasta, "ctg", _build_synthetic_genome())
    config = _make_config(tmp_path, fasta, [])
    with pytest.raises(ValueError, match=r"[Ff][Aa][Ss][Tt][Qq]"):
        run_pipeline(config)


@requires_aligners
def test_full_pipeline_recovers_known_isomer_mixture(tmp_path):
    genome = _build_synthetic_genome()
    ssc_start = LSC_LEN + IR_LEN  # 28_000
    ssc_end = ssc_start + SSC_LEN  # 34_000

    fasta = tmp_path / "genome.fasta"
    _write_fasta(fasta, "synthetic_ctg", genome)

    pairs = _simulate_reads(genome, ssc_start, ssc_end)
    fq1 = tmp_path / "reads_R1.fastq"
    fq2 = tmp_path / "reads_R2.fastq"
    _write_fastq(fq1, pairs, first=True)
    _write_fastq(fq2, pairs, first=False)

    config = _make_config(tmp_path, fasta, [fq1, fq2])
    result = run_pipeline(config)

    # ---- artifacts ---------------------------------------------------------
    assert result.bam_path is not None and result.bam_path.exists()
    assert result.html_report.exists() and result.html_report.stat().st_size > 10_000
    assert result.summary_json.exists()
    summary = json.loads(result.summary_json.read_text())
    assert summary["isomer"]["method"].startswith("dual_reference_remap")

    gfa_text = result.gfa_graph.read_text()
    assert "\tRC:i:" in gfa_text  # real segments carry assigned-read counts
    segment_lines = [line for line in gfa_text.splitlines() if line.startswith("S\t")]
    assert len(segment_lines) == 2
    assert all(len(line.split("\t")[2]) >= len(genome) * 0.9 for line in segment_lines)

    # ---- ground truth recovery --------------------------------------------
    expected_fraction_b = 300 / (700 + 300)
    isomers = result.isomer_result
    assert isomers.remap_fraction_b == pytest.approx(expected_fraction_b, abs=0.06)
    assert isomers.remap_fraction_a == pytest.approx(1 - expected_fraction_b, abs=0.06)
    # both reference frames must agree
    assert isomers.remap_agreement_delta is not None and isomers.remap_agreement_delta < 0.05
    # perfect simulated reads -> no true variants; caller must stay silent
    assert result.heteroplasmy_calls == []
    assert isomers.ssc_source == "auto_detected"
