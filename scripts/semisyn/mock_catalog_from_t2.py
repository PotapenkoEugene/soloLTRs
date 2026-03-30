"""
mock_catalog_from_t2.py — Generate a semi-synthetic catalog TSV from existing T2 ground truth.

Used to develop and test the Phase B pipeline (extract -> partition -> compose -> simulate
-> run -> evaluate) without needing to run external annotation tools (EDTA, Look4LTRs,
soloLTRseeker) on real genomes.

The T2 ground truth JSON contains synthetic solo LTRs and intact LTR-RTs with known
positions. This script converts them into the standard catalog format used by Phase B.

Usage:
    python scripts/semisyn/mock_catalog_from_t2.py \\
        --truth data/sim/t2/ground_truth/t2_ground_truth.json \\
        --regions data/pipeline/lib/regions.tsv \\
        --out data/semisyn/catalogs/t2_mock_catalog.tsv
"""

import argparse
import csv
import json
from pathlib import Path


def load_regions(path: Path) -> dict[str, dict]:
    """Load LTR/internal region lengths from regions.tsv."""
    regions = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            regions[row["family"]] = {
                "ltr_len": int(row["ltr_len"]),
                "int_len": int(row["int_len"]),
            }
    return regions


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth", required=True, type=Path,
                        help="T2 ground truth JSON")
    parser.add_argument("--regions", required=True, type=Path,
                        help="regions.tsv from sololtrs prepare")
    parser.add_argument("--out", required=True, type=Path,
                        help="Output catalog TSV")
    args = parser.parse_args()

    regions = load_regions(args.regions)
    elements = json.loads(args.truth.read_text())

    rows = []
    solo_counters: dict[str, int] = {}
    intact_counters: dict[str, int] = {}

    for elem in elements:
        etype = elem["elem_type"]
        if etype not in ("solo_LTR", "intact_LTR_RT"):
            continue

        family = elem["family"]
        chrom = elem["chrom"]
        start = elem["start"]   # 1-based inclusive (GenomeBuilder convention)
        end = elem["end"]       # 1-based inclusive
        strand = elem["strand"]
        superfamily = elem["superfamily"]

        if etype == "solo_LTR":
            solo_counters[family] = solo_counters.get(family, 0) + 1
            locus_id = f"t2_{family}_s{solo_counters[family]:04d}"
            rows.append({
                "locus_id": locus_id,
                "type": "solo_LTR",
                "family": family,
                "superfamily": superfamily,
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "ltr_5_start": ".",
                "ltr_5_end": ".",
                "ltr_3_start": ".",
                "ltr_3_end": ".",
                "int_start": ".",
                "int_end": ".",
                "n_tools": 3,
                "edta": 1,
                "look4ltrs": 1,
                "soloLTRseeker": 1,
            })

        elif etype == "intact_LTR_RT":
            if family not in regions:
                print(f"Warning: family {family} not in regions.tsv, skipping")
                continue
            ltr_len = regions[family]["ltr_len"]
            intact_counters[family] = intact_counters.get(family, 0) + 1
            locus_id = f"t2_{family}_i{intact_counters[family]:04d}"

            # Sub-element boundaries (1-based inclusive, GenomeBuilder layout):
            # [start ... ltr_5_end][ltr_5_end+1 ... int_end][int_end+1 ... end]
            ltr_5_start = start
            ltr_5_end = start + ltr_len - 1
            ltr_3_start = end - ltr_len + 1
            ltr_3_end = end
            int_start = ltr_5_end + 1
            int_end = ltr_3_start - 1

            rows.append({
                "locus_id": locus_id,
                "type": "intact_LTR_RT",
                "family": family,
                "superfamily": superfamily,
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "ltr_5_start": ltr_5_start,
                "ltr_5_end": ltr_5_end,
                "ltr_3_start": ltr_3_start,
                "ltr_3_end": ltr_3_end,
                "int_start": int_start,
                "int_end": int_end,
                "n_tools": 3,
                "edta": 1,
                "look4ltrs": 1,
                "soloLTRseeker": 1,
            })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "locus_id", "type", "family", "superfamily",
        "chrom", "start", "end", "strand",
        "ltr_5_start", "ltr_5_end", "ltr_3_start", "ltr_3_end",
        "int_start", "int_end",
        "n_tools", "edta", "look4ltrs", "soloLTRseeker",
    ]
    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    # Summary
    solos = [r for r in rows if r["type"] == "solo_LTR"]
    intacts = [r for r in rows if r["type"] == "intact_LTR_RT"]
    print(f"Written {len(rows)} loci to {args.out}")
    print(f"  Solo LTRs:    {len(solos)}")
    print(f"  Intact LTR-RTs: {len(intacts)}")
    for fam in sorted(set(r["family"] for r in rows)):
        n_s = sum(1 for r in rows if r["family"] == fam and r["type"] == "solo_LTR")
        n_i = sum(1 for r in rows if r["family"] == fam and r["type"] == "intact_LTR_RT")
        sc = f"{n_s/n_i:.2f}" if n_i > 0 else "inf"
        print(f"  {fam}: {n_s} solo, {n_i} intact, S/C={sc}")


if __name__ == "__main__":
    main()
