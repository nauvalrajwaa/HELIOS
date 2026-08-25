"""Inverted repeat / single-copy region detection for circular organelle assemblies.

Plastid (and many mitochondrial) genomes carry two large copies of an inverted
repeat (IRa / IRb). The two single-copy arcs between them are the LSC and the
SSC; reversing the SSC arc yields the alternate structural isomer. This module
locates the IR pair directly from the assembly sequence itself, so HELIOS never
has to guess boundaries (no more arbitrary ``middle 30%`` fallback).

Coordinates are 0-based half-open (start inclusive, end exclusive) on the
linear FASTA frame. Regions may wrap around the origin because plastomes are
circular; wrapped regions are reported with ``end < start`` and every consumer
in this package handles that convention.
"""

from __future__ import annotations

from dataclasses import dataclass, field

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# Typical angiosperm plastome values, used for soft validation warnings only.
EXPECTED_IR_RANGE = (10_000, 35_000)
EXPECTED_SSC_RANGE = (8_000, 40_000)


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


@dataclass(slots=True)
class RepeatRegion:
    """A repeat copy on the linear frame (may wrap the origin)."""

    name: str
    start: int  # 0-based inclusive
    end: int  # 0-based exclusive (may be < start when wrapping the origin)
    length: int


@dataclass(slots=True)
class IrDetection:
    ir_a: RepeatRegion | None = None
    ir_b: RepeatRegion | None = None
    ssc_start: int | None = None
    ssc_end: int | None = None
    lsc_start: int | None = None
    lsc_end: int | None = None
    warnings: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return (
            self.ir_a is not None
            and self.ir_b is not None
            and self.ssc_start is not None
            and self.ssc_end is not None
        )


def _complement_base(base: str) -> str:
    return base.translate(_COMP)


def _normalize_end(raw_end: int, n: int) -> int:
    """Map an exclusive end produced by modular extension onto [0, n]."""
    if raw_end <= n:
        return raw_end
    return raw_end - n


def _interval_pieces(start: int, end: int, n: int) -> list[tuple[int, int]]:
    """Split a possibly-wrapped interval into linear pieces."""
    if start <= end:
        return [(start, end)]
    return [(start, n), (0, end)]


def _overlap_length(a_start: int, a_end: int, b_start: int, b_end: int, n: int) -> int:
    """Total overlap between two possibly-wrapped intervals."""
    total = 0
    for as_, ae in _interval_pieces(a_start, a_end, n):
        for bs, be in _interval_pieces(b_start, b_end, n):
            ov = min(ae, be) - max(as_, bs)
            if ov > 0:
                total += ov
    return total


def _extend_pair(
    seq: str,
    i: int,
    j: int,
    k: int,
    n: int,
    max_span: int,
) -> tuple[int, int, int, int]:
    """Maximally extend a seeded inverted-repeat pair.

    Seed invariant: ``seq[i:i+k] == reverse_complement(seq[j:j+k])``.
    Returns ``(a_lo, a_hi_excl, b_lo, b_hi_excl)`` where the A-side interval
    equals the reverse complement of the B-side interval. Extension uses
    modular indexing so repeats spanning the assembly origin are captured.
    """

    def bases_match(a: str, b: str) -> bool:
        if a == "N" or b == "N":
            return False
        return a == _complement_base(b)

    lo_i, hi_i = i, i + k  # A-side interval [lo_i, hi_i)
    lo_j, hi_j = j, j + k  # B-side interval [lo_j, hi_j)
    total = hi_i - lo_i

    # Extend right on A-side <-> extend left on B-side.
    while total < max_span:
        a = seq[hi_i % n]
        b = seq[(lo_j - 1) % n]
        if not bases_match(a, b):
            break
        hi_i += 1
        lo_j -= 1
        total += 1

    # Extend left on A-side <-> extend right on B-side.
    while total < max_span:
        a = seq[(lo_i - 1) % n]
        b = seq[hi_j % n]
        if not bases_match(a, b):
            break
        lo_i -= 1
        hi_j += 1
        total += 1

    return lo_i, hi_i, lo_j, hi_j


def _bounding_interval(intervals: list[tuple[int, int]], n: int) -> tuple[int, int]:
    """Smallest circular arc covering every half-open interval (may wrap).

    Returns ``(start, end_excl)`` where ``end`` may exceed ``n``; callers
    derive the length as ``end - start`` and normalise coordinates with
    ``% n`` only when reporting.
    """
    coverage = bytearray(n)
    for lo, hi in intervals:
        length = (hi - lo) % n or n
        for d in range(length):
            coverage[(lo + d) % n] = 1
    if not any(coverage):
        return 0, 0
    start = next(idx for idx in range(n) if coverage[idx])
    best_gap_len = -1
    best_gap_start = 0
    gap_start: int | None = None
    pos = start
    for _ in range(n):
        if coverage[pos % n]:
            if gap_start is not None:
                gap_len = (pos % n) - gap_start
                if gap_len > best_gap_len:
                    best_gap_len, best_gap_start = gap_len, gap_start
                gap_start = None
        elif gap_start is None:
            gap_start = pos % n
        pos += 1
    if gap_start is not None:
        gap_len = start + n - gap_start
        if gap_len > best_gap_len:
            best_gap_len, best_gap_start = gap_len, gap_start
    if best_gap_len < 0:  # fully covered circle
        return 0, n
    span = n - best_gap_len
    box_start = (best_gap_start + best_gap_len) % n
    return box_start, box_start + span


def _cluster_inverted_pairs(
    pairs: list[tuple[int, int, int, int]], n: int, max_gap: int = 600
) -> list[tuple[int, int, int, int]]:
    """Chain partial repeat hits into full-length IR copies.

    Real plastome IRs often match perfectly only in patches separated by
    short diverged spacers. Two hits belong to the same copy pair when their
    A-side intervals are collinear AND their B-side intervals are collinear
    in the inverted direction, with circular gaps below ``max_gap``. Each
    cluster collapses to the bounding interval on both sides.
    """
    count = len(pairs)
    parent = list(range(count))

    def find(node: int) -> int:
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    def _pieces(iv: tuple[int, int]) -> list[tuple[int, int]]:
        """Split a possibly-wrapping half-open interval into linear pieces."""
        lo, hi = iv[0] % n, iv[1] % n
        if lo < hi:
            return [(lo, hi)]
        if lo == hi:
            return [(0, n)]
        return [(lo, n), (0, hi)]

    def _iv_dist(iv1: tuple[int, int], iv2: tuple[int, int]) -> int:
        """Circular distance between two intervals (0 when overlapping)."""
        best = n
        for s1, e1 in _pieces(iv1):
            for s2, e2 in _pieces(iv2):
                for shift in (-n, 0, n):
                    a_s, a_e = s1 + shift, e1 + shift
                    gap = max(s2 - a_e, a_s - e2)
                    best = min(best, max(gap, 0))
        return best

    def consistent(px: tuple[int, int, int, int], py: tuple[int, int, int, int]) -> bool:
        return (
            _iv_dist((px[0], px[1]), (py[0], py[1])) <= max_gap
            and _iv_dist((px[2], px[3]), (py[2], py[3])) <= max_gap
        )

    for x in range(count):
        for y in range(x + 1, count):
            if consistent(pairs[x], pairs[y]):
                rx, ry = find(x), find(y)
                if rx != ry:
                    parent[rx] = ry

    groups: dict[int, list[tuple[int, int, int, int]]] = {}
    for idx, pair in enumerate(pairs):
        groups.setdefault(find(idx), []).append(pair)

    result: list[tuple[int, int, int, int]] = []
    for members in groups.values():
        if len(members) == 1:
            result.append(members[0])
            continue
        a_lo, a_hi = _bounding_interval([(m[0] % n, m[1] % n) for m in members], n)
        b_lo, b_hi = _bounding_interval([(m[2] % n, m[3] % n) for m in members], n)
        span_a, span_b = a_hi - a_lo, b_hi - b_lo
        # sanity: both sides of a genuine IR copy must have similar span
        tolerance = 2 * max_gap + max(span_a, span_b) // 7
        if abs(span_a - span_b) <= tolerance:
            result.append((a_lo, a_hi, b_lo, b_hi))
        else:
            result.extend(members)  # suspicious cluster: keep raw hits instead
    return result


def detect_inverted_repeats(
    sequence: str,
    *,
    k: int = 31,
    seed_stride: int = 7,
    min_ir_len: int = 8_000,
    max_ir_len: int = 45_000,
) -> IrDetection:
    """Find the dominant inverted-repeat pair and derive SSC/LSC boundaries.

    The search hashes forward k-mers, looks up reverse-complement matches at a
    coarse stride, maximally extends every hit (circularity aware), merges
    redundant hits, and finally selects the two longest similar-length
    non-overlapping copies as the IRA/IRB pair.

    Returns an :class:`IrDetection`; ``ok`` is False (with explanatory
    warnings) when no credible IR pair exists.
    """
    seq = "".join(sequence.split()).upper()
    n = len(seq)
    detection = IrDetection()

    if n < 2 * min_ir_len + 2_000:
        detection.warnings.append(
            f"sequence too short ({n} bp) to host an inverted repeat pair of "
            f">{min_ir_len} bp per copy"
        )
        return detection

    # ---- hash every forward k-mer ----------------------------------------
    index: dict[str, list[int]] = {}
    for pos in range(n - k + 1):
        km = seq[pos : pos + k]
        if "N" in km:
            continue
        index.setdefault(km, []).append(pos)

    max_span = n // 2  # one IR copy cannot exceed half of a circular genome
    raw_pairs: list[tuple[int, int, int, int]] = []  # (a_lo, a_hi, b_lo, b_hi)

    for i in range(0, n - k + 1, seed_stride):
        already_found = False
        for a_lo, a_hi, _, _ in raw_pairs:
            if a_lo <= i < a_hi or (a_lo > a_hi and (i >= a_lo or i < a_hi)):
                already_found = True
                break
        if already_found:
            continue

        seed_rc = reverse_complement(seq[i : i + k])
        if "N" in seed_rc:
            continue
        partners = index.get(seed_rc)
        if not partners:
            continue

        best: tuple[int, int, int, int] | None = None
        for j in partners:
            if abs(j - i) < min(k * 4, 500):
                continue  # trivial near-self palindrome hits
            ext = _extend_pair(seq, i, j, k, n, max_span)
            if best is None or (ext[1] - ext[0]) > (best[1] - best[0]):
                best = ext
        if best is None or (best[1] - best[0]) < max(min_ir_len // 4, 500):
            continue
        # Reject degenerate self-overlapping composites: a genuine IR copy pair
        # occupies distinct loci, whereas circular wrap-around extension can
        # fold a candidate interval onto itself across the origin.
        if _overlap_length(best[0], best[1], best[2], best[3], n) > 0.25 * min(
            best[1] - best[0], best[3] - best[2]
        ):
            continue
        raw_pairs.append(best)

    # ---- merge redundant pairs (same copy discovered repeatedly) ----------
    merged: list[tuple[int, int, int, int]] = []
    for cand in sorted(raw_pairs, key=lambda p: p[1] - p[0], reverse=True):
        c_alo, c_ahi, c_blo, c_bhi = cand
        duplicate = False
        for m_alo, m_ahi, m_blo, m_bhi in merged:
            ov_a = _overlap_length(
                c_alo % n,
                _normalize_end(c_ahi, n),
                m_alo % n,
                _normalize_end(m_ahi, n),
                n,
            )
            ov_b = _overlap_length(
                c_blo % n,
                _normalize_end(c_bhi, n),
                m_blo % n,
                _normalize_end(m_bhi, n),
                n,
            )
            if ov_a > 0.5 * (c_ahi - c_alo) or ov_b > 0.5 * max(c_bhi - c_blo, 1):
                duplicate = True
                break
        if not duplicate:
            merged.append(cand)

    # ---- chain partial hits into full-length copy candidates --------------
    clustered = _cluster_inverted_pairs(merged, n)

    # ---- enumerate non-overlapping candidate copies -----------------------
    seen_boxes: set[tuple[int, int, int, int]] = set()
    pool: list[tuple[int, int, int, int]] = []
    for cand in clustered + merged:
        alen, blen = cand[1] - cand[0], cand[3] - cand[2]
        if min(alen, blen) < min_ir_len or max(alen, blen) > max_ir_len:
            continue
        norm = (cand[0] % n, cand[1] % n, cand[2] % n, cand[3] % n)
        if norm in seen_boxes:
            continue
        if any(_overlap_length(cand[0], cand[1], sel[0], sel[1], n) > 0 for sel in pool):
            continue
        seen_boxes.add(norm)
        pool.append(cand)

    if len(pool) < 2:
        detection.warnings.append(
            f"only {len(pool)} inverted-repeat arm(s) detected; cannot derive SSC boundaries"
        )
        return detection

    # ---- select the most biologically plausible pair ----------------------
    # Rank every candidate PAIR by how well its implied IR / SSC / LSC lengths
    # match typical angiosperm organelle architecture instead of blindly
    # taking the two longest hits (long spurious repeats outrank true IRs).
    def _pair_score(p: tuple[int, int, int, int], q: tuple[int, int, int, int]) -> float:
        ir_p, ir_q = p[1] - p[0], q[1] - q[0]
        arc_a = (q[0] - p[1]) % n
        arc_b = (p[0] - q[1]) % n
        ssc_len_, lsc_len_ = min(arc_a, arc_b), max(arc_a, arc_b)
        penalty = 0.0
        for length in (ir_p, ir_q):
            if not (EXPECTED_IR_RANGE[0] <= length <= EXPECTED_IR_RANGE[1]):
                penalty += 2.0
            penalty += 0.5 * abs(length - 22_500) / 22_500
        if not (EXPECTED_SSC_RANGE[0] <= ssc_len_ <= EXPECTED_SSC_RANGE[1]):
            penalty += 2.0
        penalty += 0.3 * abs(ssc_len_ - 19_000) / 19_000
        if not (40_000 <= lsc_len_ <= 110_000):
            penalty += 1.0
        return penalty

    best_pair: tuple[tuple[int, int, int, int], tuple[int, int, int, int]] | None = None
    best_score: float | None = None
    for x in range(len(pool)):
        for y in range(x + 1, len(pool)):
            score = _pair_score(pool[x], pool[y])
            improves = best_score is None or score < best_score
            if not improves and score == best_score and best_pair is not None:
                improves = (pool[x][1] - pool[x][0]) > (best_pair[0][1] - best_pair[0][0])
            if improves:
                best_score = score
                best_pair = (pool[x], pool[y])

    assert best_pair is not None  # guaranteed by len(pool) >= 2
    a1_lo, a1_hi, _, _ = best_pair[0]
    a2_lo, a2_hi, _, _ = best_pair[1]
    len1, len2 = a1_hi - a1_lo, a2_hi - a2_lo

    ratio = max(len1, len2) / max(min(len1, len2), 1)
    if ratio > 1.3:
        detection.warnings.append(
            f"IR copies differ in length by {ratio:.2f}x ({len1} vs {len2} bp); "
            "verify boundaries manually"
        )

    detection.ir_a = RepeatRegion("IRA", a1_lo % n, _normalize_end(a1_hi, n), len1)
    detection.ir_b = RepeatRegion("IRB", a2_lo % n, _normalize_end(a2_hi, n), len2)

    for length in (len1, len2):
        if not (EXPECTED_IR_RANGE[0] <= length <= EXPECTED_IR_RANGE[1]):
            detection.warnings.append(
                f"IR copy of {length} bp is outside the typical angiosperm range "
                f"{EXPECTED_IR_RANGE}"
            )

    # ---- derive the two single-copy arcs ----------------------------------
    # Arc X runs forward from the end of IR-copy-1 to the start of IR-copy-2.
    # Arc Y runs forward from the end of IR-copy-2 to the start of IR-copy-1.
    arc_x_len = (a2_lo - a1_hi) % n
    arc_y_len = (a1_lo - a2_hi) % n

    def _arc(start: int, length: int) -> tuple[int, int]:
        return (start % n, _normalize_end(start + length, n))

    if arc_x_len <= arc_y_len:
        detection.ssc_start, detection.ssc_end = _arc(a1_hi, arc_x_len)
        detection.lsc_start, detection.lsc_end = _arc(a2_hi, arc_y_len)
        ssc_len, lsc_len = arc_x_len, arc_y_len
    else:
        detection.ssc_start, detection.ssc_end = _arc(a2_hi, arc_y_len)
        detection.lsc_start, detection.lsc_end = _arc(a1_hi, arc_x_len)
        ssc_len, lsc_len = arc_y_len, arc_x_len

    if ssc_len <= 0 or lsc_len <= 0:
        detection.warnings.append(
            "degenerate single-copy arc detected; IR copies may be adjacent or duplicated"
        )
        detection.ssc_start = detection.ssc_end = None
        detection.lsc_start = detection.lsc_end = None
        return detection

    if not (EXPECTED_SSC_RANGE[0] <= ssc_len <= EXPECTED_SSC_RANGE[1]):
        detection.warnings.append(
            f"SSC of {ssc_len} bp is outside the typical angiosperm range {EXPECTED_SSC_RANGE}"
        )

    return detection


def invert_segment(sequence: str, start: int, end: int) -> str:
    """Return the sequence with the (possibly wrapping) segment reversed.

    ``start``/``end`` are 0-based half-open on the input frame; ``end < start``
    means the segment wraps the origin.
    """
    seq = "".join(sequence.split()).upper()
    n = len(seq)
    if start == end:
        raise ValueError("empty segment cannot be inverted")
    if 0 <= start < end <= n:
        return seq[:start] + reverse_complement(seq[start:end]) + seq[end:]
    if end < start:  # wraps the origin: middle = seq[start:] + seq[:end]
        head = seq[end:start]
        middle = seq[start:] + seq[:end]
        return head + reverse_complement(middle)
    raise ValueError(f"invalid segment bounds [{start}, {end}) for length {n}")
