"""
partition_samples.py — Assign loci to N virtual samples with target S/C ratios.

Each locus is assigned to exactly ONE sample (no replacement), so each sample
has a unique, non-overlapping set of real genomic sequences. Target S/C ratios
span a realistic range to probe the full operating envelope of the pipeline.

Strategy per family:
  - Sort loci randomly (seeded), then fill samples in order
  - For each sample, assign `n_intact` intact elements and `n_solo` solo LTRs
    such that n_solo / n_intact ≈ target_SC
  - Unfulfilled families (not enough loci) are reported but not included

Output manifest.json format:
  {
    "species": "...",
    "catalog": "...",
    "n_samples": N,
    "families": [...],
    "samples": [
      {
        "sample_id": "...",
        "target_SC": {"FAMILY": float, ...},
        "truth_SC":  {"FAMILY": float, ...},
        "n_solo":    {"FAMILY": int, ...},
        "n_intact":  {"FAMILY": int, ...},
        "loci": ["locus_id", ...]
      },
      ...
    ]
  }

Usage:
    python scripts/semisyn/partition_samples.py \\
        --catalog data/semisyn/catalogs/t2_mock_catalog.tsv \\
        --n-samples 8 \\
        --seed 42 \\
        --out data/semisyn/samples/t2_mock/manifest.json
"""

import argparse
import csv
import json
import math
import random
from pathlib import Path


# Default target S/C ratios — span 0 (all intact) to ~solo-only
DEFAULT_SC_TARGETS = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 20.0]

# Minimum loci required per family to include it in any sample
MIN_SOLO_PER_FAMILY = 2
MIN_INTACT_PER_FAMILY = 1

# How many intact elements to assign per sample per family (anchor)
INTACTS_PER_SAMPLE = 2


def load_catalog(path: Path) -> list[dict]:
    rows = []
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            rows.append(row)
    return rows


def partition(catalog: list[dict], n_samples: int, sc_targets: list[float],
              seed: int) -> tuple[list[dict], list[str]]:
    """
    Returns (samples list, warnings list).
    Each sample dict has: sample_id, target_SC, truth_SC, n_solo, n_intact, loci.
    """
    rng = random.Random(seed)

    # Group loci by family and type
    solos: dict[str, list[str]] = {}    # family -> [locus_id, ...]
    intacts: dict[str, list[str]] = {}  # family -> [locus_id, ...]
    for row in catalog:
        fam = row["family"]
        lid = row["locus_id"]
        if row["type"] == "solo_LTR":
            solos.setdefault(fam, []).append(lid)
        elif row["type"] == "intact_LTR_RT":
            intacts.setdefault(fam, []).append(lid)

    # Shuffle each family's loci independently
    for fam in solos:
        rng.shuffle(solos[fam])
    for fam in intacts:
        rng.shuffle(intacts[fam])

    # Determine which families have enough loci
    all_families = sorted(set(solos) | set(intacts))
    usable_families = []
    warnings = []
    for fam in all_families:
        n_s = len(solos.get(fam, []))
        n_i = len(intacts.get(fam, []))
        if n_s < MIN_SOLO_PER_FAMILY or n_i < MIN_INTACT_PER_FAMILY:
            warnings.append(
                f"Family {fam} skipped: only {n_s} solo + {n_i} intact loci "
                f"(min: {MIN_SOLO_PER_FAMILY} solo, {MIN_INTACT_PER_FAMILY} intact)"
            )
        else:
            usable_families.append(fam)

    if not usable_families:
        raise ValueError("No families have enough loci for partitioning. "
                         "Check MIN_SOLO_PER_FAMILY and MIN_INTACT_PER_FAMILY.")

    # Trim sc_targets to n_samples
    if len(sc_targets) < n_samples:
        # Repeat the pattern if fewer targets than samples
        sc_targets = (sc_targets * math.ceil(n_samples / len(sc_targets)))[:n_samples]
    else:
        sc_targets = sc_targets[:n_samples]

    # Running indices into each family's locus list
    solo_idx: dict[str, int] = {fam: 0 for fam in usable_families}
    intact_idx: dict[str, int] = {fam: 0 for fam in usable_families}

    samples = []
    for s_idx, target_sc in enumerate(sc_targets):
        sample_id = f"sample_{s_idx+1:02d}"
        assigned_loci = []
        target_sc_dict = {}
        truth_sc_dict = {}
        n_solo_dict = {}
        n_intact_dict = {}

        for fam in usable_families:
            # How many intact to use this sample?
            n_i = INTACTS_PER_SAMPLE
            # Clamp to remaining available
            available_intact = len(intacts[fam]) - intact_idx[fam]
            if available_intact < 1:
                warnings.append(
                    f"Sample {sample_id}: {fam} has no intact loci left, skipping family"
                )
                continue
            n_i = min(n_i, available_intact)

            # How many solos to achieve target_SC?
            n_s = round(target_sc * n_i)
            available_solo = len(solos.get(fam, [])) - solo_idx.get(fam, 0)
            n_s = min(n_s, available_solo)

            # Need at least 1 intact; solos can be 0
            if n_i < 1:
                continue

            # Assign loci
            fam_intact = intacts[fam][intact_idx[fam]: intact_idx[fam] + n_i]
            fam_solo = solos.get(fam, [])[solo_idx.get(fam, 0): solo_idx.get(fam, 0) + n_s]
            intact_idx[fam] += n_i
            solo_idx[fam] = solo_idx.get(fam, 0) + n_s

            assigned_loci.extend(fam_intact + fam_solo)
            truth_sc = n_s / n_i if n_i > 0 else 0.0
            target_sc_dict[fam] = target_sc
            truth_sc_dict[fam] = truth_sc
            n_solo_dict[fam] = n_s
            n_intact_dict[fam] = n_i

        samples.append({
            "sample_id": sample_id,
            "target_SC": target_sc_dict,
            "truth_SC": truth_sc_dict,
            "n_solo": n_solo_dict,
            "n_intact": n_intact_dict,
            "loci": assigned_loci,
        })

    return samples, warnings


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", required=True, type=Path)
    parser.add_argument("--n-samples", type=int, default=8)
    parser.add_argument("--sc-targets", type=str,
                        default=",".join(str(x) for x in DEFAULT_SC_TARGETS),
                        help="Comma-separated target S/C ratios")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--species", type=str, default="unknown")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    sc_targets = [float(x) for x in args.sc_targets.split(",")]
    catalog = load_catalog(args.catalog)

    samples, warnings = partition(catalog, args.n_samples, sc_targets, args.seed)

    for w in warnings:
        print(f"Warning: {w}")

    # Collect family list from samples
    all_families = sorted({fam for s in samples for fam in s["truth_SC"]})

    manifest = {
        "species": args.species,
        "catalog": str(args.catalog),
        "n_samples": len(samples),
        "families": all_families,
        "samples": samples,
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(manifest, indent=2))

    # Print summary
    print(f"\nPartitioned {len(catalog)} loci into {len(samples)} samples")
    print(f"Families: {', '.join(all_families)}")
    print(f"\n{'Sample':<12} {'Target S/C':>12}  Families x (solo/intact)")
    for s in samples:
        # Show one representative target (first family)
        if s["target_SC"]:
            fam0 = list(s["truth_SC"])[0]
            rep_target = s["target_SC"][fam0]
            rep_truth = s["truth_SC"][fam0]
            detail = "  ".join(
                f"{f}:{s['n_solo'].get(f,0)}/{s['n_intact'].get(f,0)}"
                for f in all_families if f in s["truth_SC"]
            )
            print(f"  {s['sample_id']:<10} target={rep_target:>5.1f}  truth={rep_truth:>5.2f}  {detail}")
    print(f"\nManifest written to {args.out}")


if __name__ == "__main__":
    main()
