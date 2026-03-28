"""
Build a collapsed consensus library from LTR-RT family sequences.

The collapsed consensus = 5'LTR + internal region (3'LTR removed).
This avoids the MAPQ=0 problem caused by two near-identical LTRs in the
full element: reads mapping to the LTR region would get MAPQ=0 because they
match equally well to both the 5' and 3' LTR. With only one LTR copy,
all LTR reads map unambiguously at high MAPQ.

Inputs:
  --ltr     ltr_only.fa      (individual LTR sequences per family)
  --int     internal_only.fa (internal region sequences per family)
  OR
  --lib     ltr_families.fa  (full elements; boundaries inferred from header metadata)

Output:
  collapsed_consensus.fa  — indexed for bwa-mem2
  regions.tsv             — family, ltr_len, int_len, int_start, int_end (in collapsed)
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path


@dataclass
class FamilyRegions:
    name: str
    ltr_len: int
    int_len: int
    # Coordinates in collapsed consensus (0-based, half-open)
    ltr_start: int = 0
    ltr_end: int = 0
    int_start: int = 0
    int_end: int = 0

    @property
    def collapsed_len(self) -> int:
        return self.ltr_len + self.int_len


def _parse_fasta(path: Path) -> dict[str, str]:
    """Parse FASTA → {name: seq}. Name = first word after '>' stripped of '#...'."""
    seqs: dict[str, list[str]] = {}
    current = None
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                raw = line[1:].split("#")[0].split()[0]
                # Strip common suffixes added by make_lib.py
                for suffix in ("_LTR", "_INT"):
                    raw = raw.replace(suffix, "")
                current = raw
                seqs.setdefault(current, [])
            elif current is not None:
                seqs[current].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def build_from_split(
    ltr_fa: Path,
    int_fa: Path,
    out_fa: Path,
    out_regions: Path,
) -> list[FamilyRegions]:
    """Build collapsed consensus from pre-split LTR + internal FASTAs."""
    ltr_seqs = _parse_fasta(ltr_fa)
    int_seqs = _parse_fasta(int_fa)

    shared = sorted(set(ltr_seqs) & set(int_seqs))
    if not shared:
        raise ValueError(
            f"No shared family names between {ltr_fa} and {int_fa}.\n"
            f"  LTR families: {sorted(ltr_seqs)}\n"
            f"  INT families: {sorted(int_seqs)}"
        )

    records: list[tuple[str, str]] = []
    regions: list[FamilyRegions] = []

    for name in shared:
        ltr = ltr_seqs[name]
        internal = int_seqs[name]
        collapsed = ltr + internal

        r = FamilyRegions(
            name=name,
            ltr_len=len(ltr),
            int_len=len(internal),
            ltr_start=0,
            ltr_end=len(ltr),
            int_start=len(ltr),
            int_end=len(ltr) + len(internal),
        )
        regions.append(r)
        records.append((f"{name} collapsed ltr={len(ltr)} int={len(internal)}", collapsed))

    _write_fasta(records, out_fa)
    _write_regions(regions, out_regions)

    return regions


def build_from_full(
    full_fa: Path,
    out_fa: Path,
    out_regions: Path,
) -> list[FamilyRegions]:
    """Build collapsed consensus from full LTR+internal+LTR FASTA.

    Extracts LTR length from header metadata: `ltr=NNN int=MMM`.
    Falls back to splitting the full sequence in half if metadata absent.
    """
    seqs = {}
    metas: dict[str, dict[str, int]] = {}
    current = None
    with open(full_fa) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                parts = line[1:].split()
                name = parts[0].split("#")[0]
                current = name
                seqs[name] = []
                meta = {}
                for p in parts[1:]:
                    if p.startswith("ltr="):
                        meta["ltr"] = int(p.split("=")[1])
                    elif p.startswith("internal="):
                        meta["int"] = int(p.split("=")[1])
                metas[name] = meta
            elif current is not None:
                seqs[current].append(line)

    records = []
    regions = []
    for name, seq_parts in seqs.items():
        full_seq = "".join(seq_parts)
        meta = metas.get(name, {})
        ltr_len = meta.get("ltr", 0)
        int_len = meta.get("int", 0)

        if ltr_len and int_len:
            # Use metadata
            collapsed = full_seq[:ltr_len + int_len]
        else:
            # Heuristic: full element = LTR + int + LTR → collapse to half
            # Detect by checking if first and last ~200 bp are similar (LTR symmetry)
            ltr_len = _infer_ltr_len(full_seq)
            int_len = len(full_seq) - 2 * ltr_len
            if int_len < 0:
                raise ValueError(f"Could not infer LTR length for family {name}")
            collapsed = full_seq[:ltr_len + int_len]

        r = FamilyRegions(
            name=name,
            ltr_len=ltr_len,
            int_len=int_len,
            ltr_start=0,
            ltr_end=ltr_len,
            int_start=ltr_len,
            int_end=ltr_len + int_len,
        )
        regions.append(r)
        records.append((f"{name} collapsed ltr={ltr_len} int={int_len}", collapsed))

    _write_fasta(records, out_fa)
    _write_regions(regions, out_regions)
    return regions


def _infer_ltr_len(seq: str, window: int = 200, min_identity: float = 0.80) -> int:
    """Infer LTR length by finding the longest matching prefix/suffix.

    Scans from TG at start matching CA at end.
    Returns estimated LTR length.
    """
    n = len(seq)
    # Try common LTR lengths: 200, 300, 400, 500, 600, 800, 1000, 1200, 1500, 1800
    candidates = list(range(200, min(n // 3, 2000), 50))
    best_len = candidates[0]
    best_id = 0.0
    for clen in candidates:
        prefix = seq[:clen]
        suffix = seq[n - clen:]
        matches = sum(a == b for a, b in zip(prefix, suffix))
        identity = matches / clen
        if identity > best_id:
            best_id = identity
            best_len = clen
        if identity < 0.60 and clen > 800:
            break  # identity degrades with divergence; stop searching
    return best_len


def _write_fasta(records: list[tuple[str, str]], path: Path, wrap: int = 60) -> None:
    with open(path, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            for i in range(0, len(seq), wrap):
                f.write(seq[i : i + wrap] + "\n")


def _write_regions(regions: list[FamilyRegions], path: Path) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["family", "ltr_len", "int_len", "ltr_start", "ltr_end", "int_start", "int_end"])
        for r in regions:
            w.writerow([r.name, r.ltr_len, r.int_len, r.ltr_start, r.ltr_end, r.int_start, r.int_end])


def load_regions(path: Path) -> dict[str, FamilyRegions]:
    regions = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            r = FamilyRegions(
                name=row["family"],
                ltr_len=int(row["ltr_len"]),
                int_len=int(row["int_len"]),
                ltr_start=int(row["ltr_start"]),
                ltr_end=int(row["ltr_end"]),
                int_start=int(row["int_start"]),
                int_end=int(row["int_end"]),
            )
            regions[r.name] = r
    return regions


def main():
    parser = argparse.ArgumentParser(description="Build collapsed consensus library")
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--split", nargs=2, metavar=("LTR_FA", "INT_FA"),
                     help="Use pre-split ltr_only.fa + internal_only.fa")
    src.add_argument("--lib", metavar="FULL_FA",
                     help="Use full LTR-RT library (5'LTR+internal+3'LTR)")
    parser.add_argument("--out-fa", required=True, help="Output collapsed FASTA")
    parser.add_argument("--out-regions", required=True, help="Output regions TSV")
    args = parser.parse_args()

    out_fa = Path(args.out_fa)
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    out_regions = Path(args.out_regions)

    if args.split:
        regions = build_from_split(Path(args.split[0]), Path(args.split[1]), out_fa, out_regions)
    else:
        regions = build_from_full(Path(args.lib), out_fa, out_regions)

    print(f"Built collapsed consensus for {len(regions)} families:")
    for r in regions:
        print(f"  {r.name:12s}  LTR={r.ltr_len:4d} bp  INT={r.int_len:5d} bp  total={r.collapsed_len:6d} bp")
    print(f"\nWrote: {out_fa}")
    print(f"       {out_regions}")


if __name__ == "__main__":
    main()
