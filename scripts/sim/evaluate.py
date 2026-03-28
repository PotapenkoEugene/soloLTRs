"""
Evaluate tool predictions against simulation ground truth.

Usage:
    python -m scripts.sim.evaluate \\
        --truth data/sim/t1/ground_truth/t1_ground_truth.gff3 \\
        --pred path/to/tool_predictions.gff3 \\
        --tolerance 50 \\
        --report results/eval_t1.tsv

A predicted element is a true positive if:
  - Same chromosome
  - Predicted start/end overlap ground truth start/end by ≥ tolerance bp
  - (Optionally) same strand

Metrics reported overall and stratified by:
  divergence_bin, family, tsd_status, tier
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Record:
    chrom: str
    start: int   # 1-based
    end: int     # 1-based
    strand: str
    attrs: dict


def parse_gff3(path: Path) -> list[Record]:
    records = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _, feat_type, start, end, _, strand, _, attr_str = parts[:9]
            attrs = {}
            for kv in attr_str.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attrs[k.strip()] = v.strip()
            records.append(
                Record(chrom=chrom, start=int(start), end=int(end), strand=strand, attrs=attrs)
            )
    return records


def overlaps(r1: Record, r2: Record, tolerance: int) -> bool:
    """True if r1 and r2 overlap on the same chrom and the overlap is ≥ tolerance bp."""
    if r1.chrom != r2.chrom:
        return False
    overlap = min(r1.end, r2.end) - max(r1.start, r2.start) + 1
    return overlap >= tolerance


def compute_metrics(
    truth: list[Record],
    preds: list[Record],
    tolerance: int,
    match_strand: bool = False,
) -> dict:
    """Return overall TP/FP/FN and per-stratum breakdown."""
    matched_truth = set()
    matched_pred = set()

    for i, t in enumerate(truth):
        for j, p in enumerate(preds):
            if overlaps(t, p, tolerance):
                if match_strand and t.strand != p.strand:
                    continue
                matched_truth.add(i)
                matched_pred.add(j)
                break   # one truth → one pred (greedy)

    tp = len(matched_truth)
    fn = len(truth) - tp
    fp = len(preds) - len(matched_pred)

    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    f1 = (
        2 * sensitivity * precision / (sensitivity + precision)
        if (sensitivity + precision) > 0
        else 0.0
    )

    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "n_truth": len(truth),
        "n_pred": len(preds),
        "sensitivity": sensitivity,
        "precision": precision,
        "f1": f1,
    }


def stratified_metrics(
    truth: list[Record],
    preds: list[Record],
    tolerance: int,
    strat_key: str,
) -> dict[str, dict]:
    """Compute metrics per stratum of strat_key in truth attrs."""
    strata = defaultdict(list)
    for t in truth:
        val = t.attrs.get(strat_key, "unknown")
        strata[val].append(t)
    results = {}
    for val, records in sorted(strata.items()):
        results[val] = compute_metrics(records, preds, tolerance)
    return results


def main():
    parser = argparse.ArgumentParser(description="Evaluate solo LTR detection against ground truth")
    parser.add_argument("--truth", required=True, help="Ground truth GFF3")
    parser.add_argument("--pred", required=True, help="Predictions GFF3")
    parser.add_argument("--tolerance", type=int, default=50, help="Overlap tolerance in bp (default: 50)")
    parser.add_argument("--strand", action="store_true", help="Require strand match")
    parser.add_argument("--report", default=None, help="Output TSV report path")
    args = parser.parse_args()

    truth = parse_gff3(Path(args.truth))
    preds = parse_gff3(Path(args.pred))

    # Filter truth to solo_LTRs only (background TEs excluded)
    truth_solo = [r for r in truth if r.attrs.get("ID", "").startswith("soloLTR")]
    if not truth_solo:
        truth_solo = truth  # fall back to all records

    print(f"Truth: {len(truth_solo)} solo LTRs | Predictions: {len(preds)}")
    print(f"Tolerance: ±{args.tolerance} bp\n")

    overall = compute_metrics(truth_solo, preds, args.tolerance, args.strand)
    print("Overall:")
    print(f"  TP={overall['tp']}  FP={overall['fp']}  FN={overall['fn']}")
    print(f"  Sensitivity: {overall['sensitivity']:.3f}")
    print(f"  Precision:   {overall['precision']:.3f}")
    print(f"  F1:          {overall['f1']:.3f}\n")

    rows = [{"stratum_key": "overall", "stratum_value": "all", **overall}]

    for key in ("divergence_bin", "family", "tsd_status"):
        strat = stratified_metrics(truth_solo, preds, args.tolerance, key)
        if strat:
            print(f"By {key}:")
            for val, m in strat.items():
                print(
                    f"  {val:15s}  sensitivity={m['sensitivity']:.3f}  "
                    f"precision={m['precision']:.3f}  f1={m['f1']:.3f}  "
                    f"(n={m['n_truth']})"
                )
                rows.append({"stratum_key": key, "stratum_value": val, **m})
            print()

    if args.report:
        rpath = Path(args.report)
        rpath.parent.mkdir(parents=True, exist_ok=True)
        fieldnames = list(rows[0].keys())
        with open(rpath, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        print(f"Report written to {args.report}")

    # Exit code: 0 if F1 > 0, 1 if completely empty predictions
    sys.exit(0 if overall["f1"] > 0 else 1)


if __name__ == "__main__":
    main()
