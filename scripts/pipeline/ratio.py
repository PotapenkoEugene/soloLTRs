"""
Compute solo:intact ratio (S/C) from LTR/internal depth measurements.

Math:
  R = D_INT / D_LTR = C / (S + 2C)
  S/C = (1 - 2R) / R = 1/R - 2

where S = solo LTR copies, C = intact LTR-RT copies.

R = 0.50  → S/C = 0    (all intact)
R = 0.33  → S/C = 1    (equal)
R = 0.10  → S/C = 8
R = 0.07  → S/C ≈ 12   (barley BARE-1 expected)

Binary PAV call: solo_detected = (R < R_SOLO_THRESHOLD = 0.45)
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from scripts.pipeline.config import (
    CI_ALPHA,
    MIN_LTR_DEPTH,
    MIN_MAPPED_READS,
    N_BOOTSTRAP,
    R_MAX,
    R_MIN,
    R_SOLO_THRESHOLD,
)
from scripts.pipeline.depth import DepthResult, load_depth_table


@dataclass
class RatioResult:
    family: str
    R: float                # D_INT / D_LTR
    SC: float               # solo:intact ratio (inf if C=0)
    call: str               # "all_intact" | "solo_detected" | "near_all_solo" |
                            # "low_coverage" | "absent"
    solo_detected: bool     # binary PAV
    d_ltr: float
    d_int: float
    ci_lo: float = float("nan")   # 95% CI on R (bootstrap)
    ci_hi: float = float("nan")
    sc_ci_lo: float = float("nan")
    sc_ci_hi: float = float("nan")
    confidence: str = "ok"  # "ok" | "low_reads" | "low_coverage" | "absent"


def _r_to_sc(R: float) -> float:
    """Convert depth ratio R to solo:intact ratio S/C."""
    if R <= 0:
        return float("inf")
    return (1.0 - 2.0 * R) / R


def estimate_ratio(depth: DepthResult, rng: np.random.Generator | None = None) -> RatioResult:
    """Estimate S/C from a single DepthResult."""
    d_ltr = depth.d_ltr
    d_int = depth.d_int

    # --- Guard conditions ---
    if d_ltr < 1e-6 and d_int < 1e-6:
        return RatioResult(
            family=depth.family, R=float("nan"), SC=float("nan"),
            call="absent", solo_detected=False,
            d_ltr=d_ltr, d_int=d_int, confidence="absent",
        )

    if d_ltr < MIN_LTR_DEPTH:
        return RatioResult(
            family=depth.family, R=float("nan"), SC=float("nan"),
            call="low_coverage", solo_detected=False,
            d_ltr=d_ltr, d_int=d_int, confidence="low_coverage",
        )

    R = d_int / d_ltr

    # Clamp R to [0, R_MAX]
    R_clamped = min(max(R, 0.0), R_MAX)

    if R_clamped >= R_MAX:
        call = "all_intact"
        SC = 0.0
    elif R_clamped <= R_MIN:
        call = "near_all_solo"
        SC = _r_to_sc(R_MIN)
    else:
        SC = _r_to_sc(R_clamped)
        call = "solo_detected" if R_clamped < R_SOLO_THRESHOLD else "all_intact"

    solo_detected = R_clamped < R_SOLO_THRESHOLD
    confidence = "low_reads" if depth.n_reads_ltr < MIN_MAPPED_READS else "ok"

    return RatioResult(
        family=depth.family,
        R=R_clamped,
        SC=SC,
        call=call,
        solo_detected=solo_detected,
        d_ltr=d_ltr,
        d_int=d_int,
        confidence=confidence,
    )


def bootstrap_ci(
    depth: DepthResult,
    n_boot: int = N_BOOTSTRAP,
    alpha: float = CI_ALPHA,
    seed: int = 42,
) -> tuple[float, float, float, float]:
    """Bootstrap CI on R and derived S/C.

    Treats depth values as Poisson counts: resample D_LTR and D_INT
    from Poisson(observed_depth). Returns (R_lo, R_hi, SC_lo, SC_hi).
    """
    rng = np.random.default_rng(seed)

    if depth.d_ltr < 1e-6:
        nan = float("nan")
        return nan, nan, nan, nan

    # Sample Poisson draws for both depths
    d_ltr_samples = rng.poisson(depth.d_ltr, n_boot).astype(float)
    d_int_samples = rng.poisson(depth.d_int, n_boot).astype(float)

    # Replace zeros in LTR to avoid division by zero
    d_ltr_samples = np.maximum(d_ltr_samples, 1e-6)

    R_samples = d_int_samples / d_ltr_samples
    R_samples = np.clip(R_samples, 0.0, R_MAX)

    SC_samples = np.where(
        R_samples >= R_MAX, 0.0,
        np.where(R_samples <= 0, float("inf"), (1.0 - 2.0 * R_samples) / R_samples)
    )

    lo = alpha / 2
    hi = 1.0 - lo
    R_lo, R_hi = float(np.quantile(R_samples, lo)), float(np.quantile(R_samples, hi))
    SC_hi, SC_lo = float(np.quantile(SC_samples, lo)), float(np.quantile(SC_samples, hi))
    # Note: SC is inversely related to R, so quantiles are flipped

    return R_lo, R_hi, SC_lo, SC_hi


def write_ratios(results: list[RatioResult], path: Path) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "family", "R", "SC_ratio", "call", "solo_detected",
            "D_LTR", "D_INT",
            "R_ci_lo", "R_ci_hi", "SC_ci_lo", "SC_ci_hi",
            "confidence",
        ])
        for r in results:
            def fmt(v):
                return f"{v:.4f}" if not (isinstance(v, float) and (v != v)) else "NA"
            w.writerow([
                r.family, fmt(r.R), fmt(r.SC), r.call,
                "1" if r.solo_detected else "0",
                fmt(r.d_ltr), fmt(r.d_int),
                fmt(r.ci_lo), fmt(r.ci_hi), fmt(r.sc_ci_lo), fmt(r.sc_ci_hi),
                r.confidence,
            ])


def load_ratios(path: Path) -> list[RatioResult]:
    results = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            def fv(k):
                v = row.get(k, "NA")
                return float(v) if v not in ("NA", "") else float("nan")
            results.append(RatioResult(
                family=row["family"],
                R=fv("R"), SC=fv("SC_ratio"),
                call=row["call"],
                solo_detected=row["solo_detected"] == "1",
                d_ltr=fv("D_LTR"), d_int=fv("D_INT"),
                ci_lo=fv("R_ci_lo"), ci_hi=fv("R_ci_hi"),
                sc_ci_lo=fv("SC_ci_lo"), sc_ci_hi=fv("SC_ci_hi"),
                confidence=row.get("confidence", "ok"),
            ))
    return results


def main():
    parser = argparse.ArgumentParser(description="Estimate solo:intact ratios from depth table")
    parser.add_argument("--depth", required=True, help="depth_table.tsv from depth.py")
    parser.add_argument("--out", required=True, help="Output ratios.tsv")
    parser.add_argument("--bootstrap", action="store_true", help="Compute bootstrap CI")
    parser.add_argument("--n-boot", type=int, default=N_BOOTSTRAP)
    args = parser.parse_args()

    depth_results = load_depth_table(Path(args.depth))
    ratio_results = []

    for d in depth_results:
        r = estimate_ratio(d)
        if args.bootstrap and d.d_ltr >= MIN_LTR_DEPTH:
            r.ci_lo, r.ci_hi, r.sc_ci_lo, r.sc_ci_hi = bootstrap_ci(d, n_boot=args.n_boot)
        ratio_results.append(r)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    write_ratios(ratio_results, out)

    print(f"\nSolo:intact ratio estimates ({len(ratio_results)} families):")
    print(f"{'family':12s}  {'R':>6s}  {'S/C':>8s}  {'call':18s}  {'PAV':>5s}  {'confidence'}")
    for r in ratio_results:
        R_str = f"{r.R:.3f}" if r.R == r.R else "  NA "
        SC_str = f"{r.SC:.2f}" if r.SC == r.SC and r.SC != float('inf') else (
            "  inf" if r.SC == float('inf') else "   NA"
        )
        pav = "YES" if r.solo_detected else " no"
        print(f"  {r.family:10s}  {R_str:>6s}  {SC_str:>8s}  {r.call:18s}  {pav:>5s}  {r.confidence}")

    n_detected = sum(1 for r in ratio_results if r.solo_detected)
    print(f"\n  Solo LTRs detected in {n_detected}/{len(ratio_results)} families")
    print(f"Wrote: {out}")


if __name__ == "__main__":
    main()
