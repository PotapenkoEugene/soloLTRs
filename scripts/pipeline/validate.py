"""
Compare pipeline M/U results against Cossu et al. 2017 expected S/C ratios.

Usage:
  python -m scripts.pipeline.validate \\
    --expected references/cossu/expected_ratios.tsv \\
    --results-dir data/cossu \\
    [--tolerance 0.15]

For each species in the expected file, reads data/cossu/{species}-mm2/results.tsv,
extracts the TOTAL row S_to_C, and compares against the paper's value.

Outputs:
  - Formatted comparison table to stdout
  - data/cossu/validation_summary.tsv
  - Exit 0 if all present species within tolerance, 1 otherwise
"""

import argparse
import csv
import math
import sys
from pathlib import Path


def load_expected(path: Path) -> list[dict]:
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append({
                "species": row["species"],
                "sra": row["sra"],
                "expected_M": int(row["expected_M"]),
                "expected_U": int(row["expected_U"]),
                "expected_SC": float(row["expected_SC"]),
            })
    return rows


def load_observed_sc(results_tsv: Path) -> float | None:
    """Read results.tsv and return the TOTAL row S_to_C, or None if missing/invalid."""
    if not results_tsv.exists():
        return None
    with open(results_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["family"] == "TOTAL":
                val = row["S_to_C"]
                if val in ("inf", "nan", ""):
                    return None
                return float(val)
    return None


def format_row(species, sra, expected_sc, observed_sc, tolerance) -> dict:
    if observed_sc is None:
        return {
            "species": species,
            "sra": sra,
            "expected_SC": f"{expected_sc:.3f}",
            "observed_SC": "PENDING",
            "abs_error": "—",
            "rel_error_pct": "—",
            "pass": "PENDING",
        }
    abs_err = abs(observed_sc - expected_sc)
    rel_err = abs_err / expected_sc if expected_sc != 0 else float("inf")
    passed = rel_err <= tolerance
    return {
        "species": species,
        "sra": sra,
        "expected_SC": f"{expected_sc:.3f}",
        "observed_SC": f"{observed_sc:.3f}",
        "abs_error": f"{abs_err:.3f}",
        "rel_error_pct": f"{rel_err * 100:.1f}%",
        "pass": "PASS" if passed else "FAIL",
    }


def main():
    parser = argparse.ArgumentParser(description="Validate pipeline S/C ratios vs Cossu 2017")
    parser.add_argument("--expected", required=True, type=Path, help="expected_ratios.tsv")
    parser.add_argument("--results-dir", required=True, type=Path, help="Directory containing {species}-mm2/ subdirs")
    parser.add_argument("--tolerance", type=float, default=0.15, help="Max allowed relative error (default 0.15 = 15%%)")
    args = parser.parse_args()

    expected = load_expected(args.expected)
    rows = []
    any_fail = False

    for entry in expected:
        sp = entry["species"]
        results_tsv = args.results_dir / f"{sp}-mm2" / "results.tsv"
        observed = load_observed_sc(results_tsv)
        row = format_row(sp, entry["sra"], entry["expected_SC"], observed, args.tolerance)
        rows.append(row)
        if row["pass"] == "FAIL":
            any_fail = True

    # Print formatted table
    col_widths = {
        "species": max(len(r["species"]) for r in rows),
        "sra": 12,
        "expected_SC": 11,
        "observed_SC": 11,
        "abs_error": 9,
        "rel_error_pct": 12,
        "pass": 7,
    }
    header = (
        f"{'species':<{col_widths['species']}}  "
        f"{'sra':<{col_widths['sra']}}  "
        f"{'expected_SC':>{col_widths['expected_SC']}}  "
        f"{'observed_SC':>{col_widths['observed_SC']}}  "
        f"{'abs_error':>{col_widths['abs_error']}}  "
        f"{'rel_error%':>{col_widths['rel_error_pct']}}  "
        f"{'result':<{col_widths['pass']}}"
    )
    sep = "-" * len(header)
    print(f"\nCossu et al. 2017 S/C validation (tolerance={args.tolerance*100:.0f}%)")
    print(sep)
    print(header)
    print(sep)
    for r in rows:
        print(
            f"{r['species']:<{col_widths['species']}}  "
            f"{r['sra']:<{col_widths['sra']}}  "
            f"{r['expected_SC']:>{col_widths['expected_SC']}}  "
            f"{r['observed_SC']:>{col_widths['observed_SC']}}  "
            f"{r['abs_error']:>{col_widths['abs_error']}}  "
            f"{r['rel_error_pct']:>{col_widths['rel_error_pct']}}  "
            f"{r['pass']:<{col_widths['pass']}}"
        )
    print(sep)

    # Write summary TSV
    summary_path = args.results_dir / "validation_summary.tsv"
    with open(summary_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["species", "sra", "expected_SC", "observed_SC",
                                                "abs_error", "rel_error_pct", "pass"],
                                delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nSummary written to {summary_path}")

    n_done = sum(1 for r in rows if r["pass"] != "PENDING")
    n_pass = sum(1 for r in rows if r["pass"] == "PASS")
    n_pending = sum(1 for r in rows if r["pass"] == "PENDING")
    print(f"Results: {n_pass}/{n_done} passed, {n_pending} pending\n")

    sys.exit(1 if any_fail else 0)


if __name__ == "__main__":
    main()
