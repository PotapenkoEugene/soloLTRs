"""
Core utilities for building synthetic solo LTR simulation datasets.

Classes / functions:
  generate_background(length, gc_content, seed)  → str
  mutate_sequence(seq, divergence, seed)          → str
  make_tsd(length, seed)                          → str
  degrade_tsd(tsd, mode, seed)                    → str
  GenomeBuilder                                   — coordinate-tracking genome builder
  write_gff3(elements, path)
  write_bed(elements, path)
  write_manifest(elements, path)                  — JSON per-element metadata
"""

from __future__ import annotations

import json
import random
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Literal


# ---------------------------------------------------------------------------
# Sequence generation helpers
# ---------------------------------------------------------------------------

def generate_background(length: int, gc_content: float = 0.42, seed: int = 0) -> str:
    """Return a random DNA sequence with controlled GC content."""
    rng = random.Random(seed)
    at = 1.0 - gc_content
    pool = (
        "G" * int(gc_content / 2 * 200)
        + "C" * int(gc_content / 2 * 200)
        + "A" * int(at / 2 * 200)
        + "T" * int(at / 2 * 200)
    )
    return "".join(rng.choices(list(pool), k=length))


def mutate_sequence(seq: str, divergence: float, seed: int = 0) -> str:
    """Apply random substitutions at given rate (Jukes-Cantor model).

    divergence: fraction of sites to mutate (0.0–1.0).
    Each mutated site picks uniformly from the 3 other nucleotides.
    """
    if divergence <= 0:
        return seq
    rng = random.Random(seed)
    bases = list("ACGT")
    seq_list = list(seq)
    n_subs = max(1, int(len(seq) * divergence))
    positions = rng.sample(range(len(seq)), k=min(n_subs, len(seq)))
    for pos in positions:
        orig = seq_list[pos].upper()
        alts = [b for b in bases if b != orig]
        if alts:
            seq_list[pos] = rng.choice(alts)
    return "".join(seq_list)


def make_tsd(length: int = 5, seed: int = 0) -> str:
    """Generate a random TSD sequence."""
    rng = random.Random(seed)
    return "".join(rng.choices("ACGT", k=length))


def degrade_tsd(
    tsd: str,
    mode: Literal["perfect", "mismatch_1", "degraded"] = "perfect",
    seed: int = 0,
) -> tuple[str, str]:
    """Return (tsd_5, tsd_3) with controlled degradation.

    perfect     → tsd_5 == tsd_3
    mismatch_1  → one position differs
    degraded    → 1-2 bp deletion in tsd_3
    """
    if mode == "perfect":
        return tsd, tsd

    rng = random.Random(seed)

    if mode == "mismatch_1":
        pos = rng.randrange(len(tsd))
        alts = [b for b in "ACGT" if b != tsd[pos]]
        tsd_3 = tsd[:pos] + rng.choice(alts) + tsd[pos + 1 :]
        return tsd, tsd_3

    # degraded: 1-2 bp deletion in tsd_3
    n_del = rng.randint(1, min(2, len(tsd) - 1))
    del_start = rng.randint(0, len(tsd) - n_del)
    tsd_3 = tsd[:del_start] + tsd[del_start + n_del :]
    return tsd, tsd_3


# ---------------------------------------------------------------------------
# Element record
# ---------------------------------------------------------------------------

@dataclass
class Element:
    """Represents a single TE element placed in the genome."""
    elem_id: str
    elem_type: Literal["solo_LTR", "intact_LTR_RT", "background_TE"]
    family: str
    superfamily: str
    chrom: str
    start: int           # 1-based, inclusive
    end: int             # 1-based, inclusive
    strand: str          # "+" or "-"
    # Structural metadata
    divergence: float = 0.0
    divergence_bin: str = ""   # "0-5", "5-10", "10-20", "20-30"
    tsd_5: str = ""
    tsd_3: str = ""
    tsd_status: str = "perfect"   # perfect | mismatch_1 | degraded
    # Context metadata
    tier: int = 0
    nested_in: str = ""       # parent element ID if nested
    cluster_id: str = ""      # tandem cluster ID if applicable
    is_truncated_confounder: bool = False
    is_gce_pair: str = ""     # partner element ID for GCE pairs


# ---------------------------------------------------------------------------
# GenomeBuilder
# ---------------------------------------------------------------------------

class GenomeBuilder:
    """Build a synthetic genome by inserting elements into a background sequence.

    Tracks coordinate shifts as insertions lengthen the sequence.
    Supports nested insertions (insert into an existing element).

    Usage:
        builder = GenomeBuilder(background_seq, chrom_name, tier=1)
        elem = builder.insert_solo_ltr(position, ltr_seq, family, ...)
        fasta, elements = builder.build()
    """

    def __init__(self, background: str, chrom: str = "chr1", tier: int = 0):
        self._seq = list(background)
        self._offset = 0        # cumulative shift from all insertions so far
        self._elements: list[Element] = []
        self._chrom = chrom
        self._tier = tier
        self._counter = 0

    def _next_id(self, prefix: str) -> str:
        self._counter += 1
        return f"{prefix}_{self._counter:04d}"

    def _adjusted_position(self, pos: int) -> int:
        """Translate a pre-insertion position to current sequence position."""
        return pos + self._offset

    def _do_insert(self, pos: int, seq_to_insert: str, tsd: str = "") -> tuple[int, int]:
        """Insert seq_to_insert at adjusted position pos.

        If tsd is provided, the TSD target site is duplicated on both sides.
        Returns (start_1based, end_1based) in final coordinates.
        """
        adj = pos
        if tsd:
            # Duplicate tsd_5 on the left, insert element, tsd_3 already matched
            insert = tsd + seq_to_insert + tsd
        else:
            insert = seq_to_insert
        self._seq[adj:adj] = list(insert)
        start = adj + (len(tsd) if tsd else 0) + 1   # 1-based
        end = start + len(seq_to_insert) - 1
        self._offset += len(insert)
        return start, end

    def insert_solo_ltr(
        self,
        position: int,
        ltr_seq: str,
        family: str,
        superfamily: str,
        strand: str = "+",
        divergence: float = 0.0,
        tsd: str = "",
        tsd_status: str = "perfect",
        tsd_3: str = "",
        divergence_bin: str = "",
        seed: int = 0,
        tier: int | None = None,
    ) -> Element:
        """Insert a solo LTR at the given (pre-shift) position."""
        if strand == "-":
            ltr_seq = _revcomp(ltr_seq)
        start, end = self._do_insert(position + self._offset, ltr_seq, tsd)
        elem = Element(
            elem_id=self._next_id("soloLTR"),
            elem_type="solo_LTR",
            family=family,
            superfamily=superfamily,
            chrom=self._chrom,
            start=start,
            end=end,
            strand=strand,
            divergence=divergence,
            divergence_bin=divergence_bin or _div_bin(divergence),
            tsd_5=tsd,
            tsd_3=tsd_3 if tsd_3 else tsd,
            tsd_status=tsd_status,
            tier=tier if tier is not None else self._tier,
        )
        self._elements.append(elem)
        return elem

    def insert_intact_ltr_rt(
        self,
        position: int,
        ltr_seq: str,
        internal_seq: str,
        family: str,
        superfamily: str,
        strand: str = "+",
        divergence: float = 0.0,
        tsd: str = "",
        tsd_status: str = "perfect",
        tsd_3: str = "",
    ) -> Element:
        """Insert a full LTR-RT: 5'LTR + internal + 3'LTR."""
        ltr_5 = ltr_seq
        ltr_3 = mutate_sequence(ltr_seq, divergence, seed=hash(ltr_seq[:10]))
        internal = internal_seq
        full = ltr_5 + internal + ltr_3
        if strand == "-":
            full = _revcomp(full)
        start, end = self._do_insert(position + self._offset, full, tsd)
        elem = Element(
            elem_id=self._next_id("intact"),
            elem_type="intact_LTR_RT",
            family=family,
            superfamily=superfamily,
            chrom=self._chrom,
            start=start,
            end=end,
            strand=strand,
            divergence=divergence,
            divergence_bin=_div_bin(divergence),
            tsd_5=tsd,
            tsd_3=tsd_3 if tsd_3 else tsd,
            tsd_status=tsd_status,
            tier=self._tier,
        )
        self._elements.append(elem)
        return elem

    def insert_background_te(
        self,
        position: int,
        te_seq: str,
        family: str,
        strand: str = "+",
    ) -> Element:
        """Insert a non-LTR background TE fragment."""
        if strand == "-":
            te_seq = _revcomp(te_seq)
        start, end = self._do_insert(position + self._offset, te_seq)
        elem = Element(
            elem_id=self._next_id("bgTE"),
            elem_type="background_TE",
            family=family,
            superfamily="non-LTR",
            chrom=self._chrom,
            start=start,
            end=end,
            strand=strand,
            tier=self._tier,
        )
        self._elements.append(elem)
        return elem

    def get_sequence(self) -> str:
        return "".join(self._seq)

    def get_elements(self) -> list[Element]:
        return list(self._elements)


# ---------------------------------------------------------------------------
# Ground truth writers
# ---------------------------------------------------------------------------

def write_fasta(records: list[tuple[str, str]], path: Path, wrap: int = 60) -> None:
    """Write FASTA file. records = [(header, seq), ...]"""
    with open(path, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            for i in range(0, len(seq), wrap):
                f.write(seq[i : i + wrap] + "\n")


def write_gff3(elements: list[Element], path: Path) -> None:
    """Write GFF3 ground truth file."""
    with open(path, "w") as f:
        f.write("##gff-version 3\n")
        for e in elements:
            if e.elem_type not in ("solo_LTR", "intact_LTR_RT"):
                continue
            feature_type = "solo_LTR" if e.elem_type == "solo_LTR" else "LTR_retrotransposon"
            attrs = _gff3_attrs(e)
            f.write(
                f"{e.chrom}\tsim\t{feature_type}\t{e.start}\t{e.end}\t.\t"
                f"{e.strand}\t.\t{attrs}\n"
            )


def write_bed(elements: list[Element], path: Path) -> None:
    """Write BED6 ground truth file (0-based half-open coordinates)."""
    with open(path, "w") as f:
        for e in elements:
            if e.elem_type not in ("solo_LTR", "intact_LTR_RT"):
                continue
            start0 = e.start - 1
            end0 = e.end
            f.write(
                f"{e.chrom}\t{start0}\t{end0}\t{e.elem_id}\t0\t{e.strand}\n"
            )


def write_manifest(elements: list[Element], path: Path) -> None:
    """Write JSON manifest with full per-element metadata."""
    records = [asdict(e) for e in elements]
    with open(path, "w") as f:
        json.dump(records, f, indent=2)


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def _div_bin(div: float) -> str:
    if div < 0.05:
        return "0-5"
    if div < 0.10:
        return "5-10"
    if div < 0.20:
        return "10-20"
    return "20-30"


def _gff3_attrs(e: Element) -> str:
    parts = [
        f"ID={e.elem_id}",
        f"family={e.family}",
        f"superfamily={e.superfamily}",
        f"divergence={e.divergence:.4f}",
        f"divergence_bin={e.divergence_bin}",
        f"tsd_5={e.tsd_5 or 'NA'}",
        f"tsd_3={e.tsd_3 or 'NA'}",
        f"tsd_status={e.tsd_status}",
        f"tier={e.tier}",
    ]
    if e.nested_in:
        parts.append(f"nested_in={e.nested_in}")
    if e.cluster_id:
        parts.append(f"cluster_id={e.cluster_id}")
    return ";".join(parts)
