"""
Aggregate per-sample ratio results into population-level matrices.

Inputs:  directory of per-sample ratios.tsv files (or a samples.tsv manifest)
Outputs:
  ratio_matrix.tsv   — families × samples (S/C continuous values)
  R_matrix.tsv       — families × samples (R = D_INT/D_LTR continuous)
  pav_matrix.tsv     — families × samples (binary: 1=solo_detected, 0=not)
  summary.tsv        — per-family population statistics
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.pipeline.ratio import load_ratios


def build_matrices(
    sample_ratios: dict[str, Path],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return (SC_matrix, R_matrix, PAV_matrix) as DataFrames (families × samples)."""
    all_data: dict[str, dict[str, dict]] = {}  # sample → family → RatioResult

    families: set[str] = set()
    for sample_id, path in sample_ratios.items():
        results = load_ratios(path)
        all_data[sample_id] = {r.family: r for r in results}
        families.update(r.family for r in results)

    families = sorted(families)
    samples = sorted(sample_ratios)

    sc_matrix = pd.DataFrame(index=families, columns=samples, dtype=float)
    r_matrix = pd.DataFrame(index=families, columns=samples, dtype=float)
    pav_matrix = pd.DataFrame(index=families, columns=samples, dtype=int)

    for sample in samples:
        for fam in families:
            result = all_data[sample].get(fam)
            if result is None or result.confidence in ("absent", "low_coverage"):
                sc_matrix.loc[fam, sample] = float("nan")
                r_matrix.loc[fam, sample] = float("nan")
                pav_matrix.loc[fam, sample] = -1  # missing
            else:
                sc_matrix.loc[fam, sample] = result.SC
                r_matrix.loc[fam, sample] = result.R
                pav_matrix.loc[fam, sample] = 1 if result.solo_detected else 0

    return sc_matrix, r_matrix, pav_matrix


def build_summary(sc_matrix: pd.DataFrame, pav_matrix: pd.DataFrame) -> pd.DataFrame:
    """Per-family population statistics."""
    rows = []
    for fam in sc_matrix.index:
        sc_vals = sc_matrix.loc[fam].dropna()
        pav_vals = pav_matrix.loc[fam]
        n_assessed = int((pav_vals >= 0).sum())
        n_detected = int((pav_vals == 1).sum())
        rows.append({
            "family": fam,
            "n_samples_assessed": n_assessed,
            "n_solo_detected": n_detected,
            "pct_solo_detected": round(100 * n_detected / n_assessed, 1) if n_assessed > 0 else float("nan"),
            "mean_SC": round(float(sc_vals.mean()), 3) if len(sc_vals) > 0 else float("nan"),
            "median_SC": round(float(sc_vals.median()), 3) if len(sc_vals) > 0 else float("nan"),
            "std_SC": round(float(sc_vals.std()), 3) if len(sc_vals) > 1 else float("nan"),
            "min_SC": round(float(sc_vals.min()), 3) if len(sc_vals) > 0 else float("nan"),
            "max_SC": round(float(sc_vals.max()), 3) if len(sc_vals) > 0 else float("nan"),
        })
    return pd.DataFrame(rows).set_index("family")


def write_matrix(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", na_rep="NA")


def main():
    parser = argparse.ArgumentParser(description="Build population PAV matrix from per-sample ratios")
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--results-dir", metavar="DIR",
                     help="Directory containing per-sample subdirs with ratios.tsv")
    src.add_argument("--samples", metavar="TSV",
                     help="TSV with columns: sample_id, ratios_tsv")
    parser.add_argument("--out", required=True, help="Output directory")
    args = parser.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    sample_ratios: dict[str, Path] = {}

    if args.results_dir:
        results_dir = Path(args.results_dir)
        for subdir in sorted(results_dir.iterdir()):
            if subdir.is_dir():
                ratios_file = subdir / "ratios.tsv"
                if ratios_file.exists():
                    sample_ratios[subdir.name] = ratios_file
    else:
        with open(args.samples) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                sample_ratios[row["sample_id"]] = Path(row["ratios_tsv"])

    if not sample_ratios:
        import sys
        sys.exit("ERROR: No sample ratios found.")

    print(f"Building matrices for {len(sample_ratios)} samples...")
    sc_matrix, r_matrix, pav_matrix = build_matrices(sample_ratios)
    summary = build_summary(sc_matrix, pav_matrix)

    write_matrix(sc_matrix, out / "SC_matrix.tsv")
    write_matrix(r_matrix, out / "R_matrix.tsv")
    write_matrix(pav_matrix, out / "pav_matrix.tsv")
    write_matrix(summary, out / "summary.tsv")

    print(f"\nPopulation summary ({len(sc_matrix.index)} families × {len(sc_matrix.columns)} samples):")
    print(summary.to_string())
    print(f"\nWrote: {out}/")


if __name__ == "__main__":
    main()
