"""
Compute per-region mean depth for each LTR-RT family from a BAM file.

For each family in regions.tsv:
  - D_LTR = mean depth across the LTR region (trimmed)
  - D_INT = mean depth across the internal region (trimmed)

Edge trimming (default 50 bp) removes boundary artifacts where reads
spanning the LTR→internal junction have split alignments.

Depth is computed by parsing `samtools depth` output directly in Python
(no pysam dependency required).
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from scripts.pipeline.annotate import FamilyRegions, load_regions
from scripts.pipeline.config import EDGE_TRIM_BP, MIN_MAPQ_DEPTH


@dataclass
class DepthResult:
    family: str
    d_ltr: float          # mean depth in LTR core (trimmed)
    d_int: float          # mean depth in internal core (trimmed)
    ltr_len: int          # trimmed LTR length used
    int_len: int          # trimmed internal length used
    n_reads_ltr: int      # reads covering LTR region
    n_reads_int: int      # reads covering internal region
    ltr_covered_frac: float  # fraction of LTR positions with >0 depth
    int_covered_frac: float


def _samtools_depth(
    bam: Path,
    ref_name: str,
    start: int,
    end: int,
    min_mapq: int = MIN_MAPQ_DEPTH,
) -> np.ndarray:
    """Run samtools depth over [start, end) on ref_name. Returns per-position array."""
    # samtools depth uses 1-based coordinates
    cmd = [
        "samtools", "depth",
        "-a",                    # output all positions (including zero depth)
        "-q", str(min_mapq),     # minimum base quality
        "-Q", str(min_mapq),     # minimum mapping quality
        "-r", f"{ref_name}:{start + 1}-{end}",
        str(bam),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f"ERROR: samtools depth failed:\n{result.stderr}")

    depths = []
    for line in result.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 3:
            depths.append(int(parts[2]))

    if not depths:
        return np.zeros(end - start, dtype=np.float32)

    arr = np.array(depths, dtype=np.float32)
    # Pad/trim to exact length
    expected_len = end - start
    if len(arr) < expected_len:
        arr = np.pad(arr, (0, expected_len - len(arr)))
    elif len(arr) > expected_len:
        arr = arr[:expected_len]
    return arr


def compute_depth(
    bam: Path,
    regions: dict[str, FamilyRegions],
    trim: int = EDGE_TRIM_BP,
) -> list[DepthResult]:
    """Compute depth for all families. trim = bp to exclude at each region edge."""
    results = []

    for name, r in regions.items():
        # Trimmed LTR region: [ltr_start + trim, ltr_end - trim)
        ltr_core_start = r.ltr_start + trim
        ltr_core_end = r.ltr_end - trim
        if ltr_core_end <= ltr_core_start:
            # LTR too short to trim — use full region
            ltr_core_start = r.ltr_start
            ltr_core_end = r.ltr_end

        # Trimmed internal region: [int_start + trim, int_end - trim)
        int_core_start = r.int_start + trim
        int_core_end = r.int_end - trim
        if int_core_end <= int_core_start:
            int_core_start = r.int_start
            int_core_end = r.int_end

        ltr_depths = _samtools_depth(bam, name, ltr_core_start, ltr_core_end)
        int_depths = _samtools_depth(bam, name, int_core_start, int_core_end)

        d_ltr = float(np.mean(ltr_depths))
        d_int = float(np.mean(int_depths))

        # Count reads: positions with depth > 0
        ltr_covered = float(np.mean(ltr_depths > 0))
        int_covered = float(np.mean(int_depths > 0))

        # Approximate read counts (not exact — just for reporting)
        n_reads_ltr = int(np.sum(ltr_depths) // 150)  # assuming 150 bp reads
        n_reads_int = int(np.sum(int_depths) // 150)

        results.append(DepthResult(
            family=name,
            d_ltr=d_ltr,
            d_int=d_int,
            ltr_len=ltr_core_end - ltr_core_start,
            int_len=int_core_end - int_core_start,
            n_reads_ltr=n_reads_ltr,
            n_reads_int=n_reads_int,
            ltr_covered_frac=ltr_covered,
            int_covered_frac=int_covered,
        ))

    return results


def write_depth_table(results: list[DepthResult], path: Path) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "family", "d_ltr", "d_int",
            "ltr_len", "int_len",
            "n_reads_ltr", "n_reads_int",
            "ltr_covered_frac", "int_covered_frac",
        ])
        for r in results:
            w.writerow([
                r.family,
                f"{r.d_ltr:.4f}", f"{r.d_int:.4f}",
                r.ltr_len, r.int_len,
                r.n_reads_ltr, r.n_reads_int,
                f"{r.ltr_covered_frac:.4f}", f"{r.int_covered_frac:.4f}",
            ])


def load_depth_table(path: Path) -> list[DepthResult]:
    results = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            results.append(DepthResult(
                family=row["family"],
                d_ltr=float(row["d_ltr"]),
                d_int=float(row["d_int"]),
                ltr_len=int(row["ltr_len"]),
                int_len=int(row["int_len"]),
                n_reads_ltr=int(row["n_reads_ltr"]),
                n_reads_int=int(row["n_reads_int"]),
                ltr_covered_frac=float(row["ltr_covered_frac"]),
                int_covered_frac=float(row["int_covered_frac"]),
            ))
    return results


def main():
    parser = argparse.ArgumentParser(description="Compute per-family LTR/internal depths")
    parser.add_argument("--bam", required=True, help="Sorted, indexed BAM")
    parser.add_argument("--regions", required=True, help="regions.tsv from annotate.py")
    parser.add_argument("--out", required=True, help="Output depth_table.tsv")
    parser.add_argument("--trim", type=int, default=EDGE_TRIM_BP, help="Edge trim bp")
    args = parser.parse_args()

    regions = load_regions(Path(args.regions))
    results = compute_depth(Path(args.bam), regions, trim=args.trim)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    write_depth_table(results, out)

    print(f"\nDepth summary ({len(results)} families):")
    print(f"{'family':12s}  {'D_LTR':>8s}  {'D_INT':>8s}  {'coverage':>10s}")
    for r in results:
        print(
            f"  {r.family:10s}  {r.d_ltr:8.1f}  {r.d_int:8.1f}  "
            f"LTR={r.ltr_covered_frac:.1%} INT={r.int_covered_frac:.1%}"
        )
    print(f"\nWrote: {out}")


if __name__ == "__main__":
    main()
