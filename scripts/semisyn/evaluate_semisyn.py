"""
evaluate_semisyn.py — Compare pipeline S/C estimates against semi-synthetic ground truth.

For each sample x coverage combination, reads the pipeline's ratios.tsv and the
per-sample truth.json, then computes per-family accuracy metrics.

Metrics per (species, sample, coverage, family):
  - SC_truth         : true solo:intact ratio
  - SC_estimated     : pipeline estimate
  - error            : SC_estimated - SC_truth
  - abs_error        : |error|
  - relative_error   : error / (SC_truth + 1)  (avoids division-by-zero)
  - within_CI        : SC_truth in [SC_ci_lo, SC_ci_hi]
  - call             : pipeline call (solo_detected, all_intact, low_coverage, ...)

Aggregate metrics (eval_summary.tsv):
  - RMSE, MAE, bias, Pearson r, Spearman rho, CI calibration (% truth in 95% CI)
  Stratified by: species, coverage, family

Usage:
    python scripts/semisyn/evaluate_semisyn.py \\
        --manifest data/semisyn/samples/t2_mock/manifest.json \\
        --results-dir data/semisyn/results/t2_mock \\
        --coverages 5,10,30 \\
        --out-dir data/semisyn/eval/t2_mock
"""

import argparse
import csv
import json
import math
import sys
from pathlib import Path


def load_ratios(path: Path) -> dict[str, dict]:
    """Load ratios.tsv -> {family: row_dict}."""
    result = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            result[row["family"]] = row
    return result


def load_truth(truth_path: Path) -> dict:
    return json.loads(truth_path.read_text())


def pearson(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    if n < 2:
        return float("nan")
    mx, my = sum(xs) / n, sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    denom = math.sqrt(sum((x - mx)**2 for x in xs) * sum((y - my)**2 for y in ys))
    return num / denom if denom > 0 else float("nan")


def spearman(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    if n < 2:
        return float("nan")
    rx = [sorted(xs).index(x) for x in xs]
    ry = [sorted(ys).index(y) for y in ys]
    return pearson(rx, ry)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--results-dir", required=True, type=Path,
                        help="Directory with sample_XX_YYx/ subdirs")
    parser.add_argument("--coverages", type=str, default="5,10,30",
                        help="Comma-separated coverage levels to look for")
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()

    coverages = [int(c) for c in args.coverages.split(",")]
    manifest = json.loads(args.manifest.read_text())
    species = manifest.get("species", "unknown")

    args.out_dir.mkdir(parents=True, exist_ok=True)

    detail_rows = []
    missing_results = []

    for sample in manifest["samples"]:
        sample_id = sample["sample_id"]
        truth_sc = sample["truth_SC"]

        if not truth_sc:
            continue  # empty sample (no loci assigned)

        for cov in coverages:
            result_dir = args.results_dir / f"{sample_id}_{cov}x"
            ratios_path = result_dir / "ratios.tsv"

            if not ratios_path.exists():
                missing_results.append(str(ratios_path))
                continue

            ratios = load_ratios(ratios_path)

            for family, sc_true in truth_sc.items():
                row = ratios.get(family)
                if row is None:
                    # Family absent from pipeline output
                    detail_rows.append({
                        "species": species,
                        "sample_id": sample_id,
                        "coverage": cov,
                        "family": family,
                        "SC_truth": sc_true,
                        "n_solo_truth": sample["n_solo"].get(family, 0),
                        "n_intact_truth": sample["n_intact"].get(family, 0),
                        "SC_estimated": "NA",
                        "error": "NA",
                        "abs_error": "NA",
                        "relative_error": "NA",
                        "within_CI": "NA",
                        "call": "absent_from_output",
                        "D_LTR": "NA",
                        "D_INT": "NA",
                        "SC_ci_lo": "NA",
                        "SC_ci_hi": "NA",
                    })
                    continue

                try:
                    sc_est = float(row["SC_ratio"])
                    ci_lo = float(row["SC_ci_lo"])
                    ci_hi = float(row["SC_ci_hi"])
                    d_ltr = float(row["D_LTR"])
                    d_int = float(row["D_INT"])
                except (ValueError, KeyError):
                    sc_est = float("nan")
                    ci_lo = ci_hi = d_ltr = d_int = float("nan")

                error = sc_est - sc_true
                abs_error = abs(error)
                rel_error = error / (sc_true + 1.0)
                within_ci = "NA"
                if not math.isnan(ci_lo) and not math.isnan(ci_hi):
                    within_ci = "1" if ci_lo <= sc_true <= ci_hi else "0"

                detail_rows.append({
                    "species": species,
                    "sample_id": sample_id,
                    "coverage": cov,
                    "family": family,
                    "SC_truth": sc_true,
                    "n_solo_truth": sample["n_solo"].get(family, 0),
                    "n_intact_truth": sample["n_intact"].get(family, 0),
                    "SC_estimated": f"{sc_est:.4f}" if not math.isnan(sc_est) else "NA",
                    "error": f"{error:.4f}" if not math.isnan(error) else "NA",
                    "abs_error": f"{abs_error:.4f}" if not math.isnan(abs_error) else "NA",
                    "relative_error": f"{rel_error:.4f}" if not math.isnan(rel_error) else "NA",
                    "within_CI": within_ci,
                    "call": row.get("call", "NA"),
                    "D_LTR": f"{d_ltr:.2f}" if not math.isnan(d_ltr) else "NA",
                    "D_INT": f"{d_int:.2f}" if not math.isnan(d_int) else "NA",
                    "SC_ci_lo": f"{ci_lo:.4f}" if not math.isnan(ci_lo) else "NA",
                    "SC_ci_hi": f"{ci_hi:.4f}" if not math.isnan(ci_hi) else "NA",
                })

    if missing_results:
        print(f"Warning: {len(missing_results)} result files not found "
              f"(run semisyn-run first)", file=sys.stderr)
        for p in missing_results[:10]:
            print(f"  {p}", file=sys.stderr)

    # Write detail TSV
    if detail_rows:
        detail_path = args.out_dir / "eval_detail.tsv"
        fieldnames = list(detail_rows[0].keys())
        with open(detail_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(detail_rows)
        print(f"Written {len(detail_rows)} rows to {detail_path}")

    # Compute summary metrics per (species, coverage)
    def compute_agg(rows: list[dict]) -> dict:
        numeric = [(float(r["SC_truth"]), float(r["SC_estimated"]))
                   for r in rows
                   if r["SC_estimated"] not in ("NA", "") and r["SC_truth"] not in ("NA", "")]
        if not numeric:
            return {"n": 0, "RMSE": "NA", "MAE": "NA", "bias": "NA",
                    "pearson_r": "NA", "spearman_rho": "NA", "CI_calibration": "NA"}
        truths = [x[0] for x in numeric]
        ests = [x[1] for x in numeric]
        errors = [e - t for t, e in numeric]
        rmse = math.sqrt(sum(e**2 for e in errors) / len(errors))
        mae = sum(abs(e) for e in errors) / len(errors)
        bias = sum(errors) / len(errors)
        ci_rows = [r for r in rows if r["within_CI"] in ("0", "1")]
        ci_cal = (sum(1 for r in ci_rows if r["within_CI"] == "1") / len(ci_rows)
                  if ci_rows else float("nan"))
        return {
            "n": len(numeric),
            "RMSE": f"{rmse:.4f}",
            "MAE": f"{mae:.4f}",
            "bias": f"{bias:.4f}",
            "pearson_r": f"{pearson(truths, ests):.4f}",
            "spearman_rho": f"{spearman(truths, ests):.4f}",
            "CI_calibration": f"{ci_cal:.3f}" if not math.isnan(ci_cal) else "NA",
        }

    summary_rows = []
    # Overall
    agg = compute_agg(detail_rows)
    summary_rows.append({"strata": "all", "species": species, "coverage": "all",
                         "family": "all", **agg})
    # Per coverage
    for cov in coverages:
        cov_rows = [r for r in detail_rows if str(r["coverage"]) == str(cov)]
        agg = compute_agg(cov_rows)
        summary_rows.append({"strata": "coverage", "species": species,
                              "coverage": cov, "family": "all", **agg})
    # Per family
    families = sorted(set(r["family"] for r in detail_rows))
    for fam in families:
        fam_rows = [r for r in detail_rows if r["family"] == fam]
        agg = compute_agg(fam_rows)
        summary_rows.append({"strata": "family", "species": species,
                              "coverage": "all", "family": fam, **agg})

    summary_path = args.out_dir / "eval_summary.tsv"
    if summary_rows:
        with open(summary_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()),
                                    delimiter="\t")
            writer.writeheader()
            writer.writerows(summary_rows)
        print(f"Written {len(summary_rows)} summary rows to {summary_path}")

    # Print key metrics to stdout
    overall = next((r for r in summary_rows if r["strata"] == "all"), None)
    if overall:
        print(f"\nOverall (n={overall['n']}): "
              f"RMSE={overall['RMSE']}  MAE={overall['MAE']}  "
              f"bias={overall['bias']}  r={overall['pearson_r']}  "
              f"CI_cal={overall['CI_calibration']}")
    for cov in coverages:
        row = next((r for r in summary_rows
                    if r["strata"] == "coverage" and str(r["coverage"]) == str(cov)), None)
        if row and row["n"] > 0:
            print(f"  {cov}x (n={row['n']}): "
                  f"RMSE={row['RMSE']}  r={row['pearson_r']}  CI_cal={row['CI_calibration']}")


if __name__ == "__main__":
    main()
