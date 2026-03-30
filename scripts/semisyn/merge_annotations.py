"""
merge_annotations.py — 3-tool consensus intersection for solo LTR / intact LTR-RT loci.

Takes per-tool loci TSVs (from parse_edta.py, parse_look4ltrs.py, parse_soloLTRseeker.py)
and produces a high-confidence catalog by requiring agreement across all three tools.

Merge logic:
  1. Cluster loci from all tools by reciprocal overlap >= MIN_OVERLAP (default 50%).
     Two loci overlap reciprocally if:
       overlap_len / min(len_a, len_b) >= MIN_OVERLAP
  2. Retain clusters where ALL 3 tools contributed at least one locus of the same type.
     - Solo clusters: all 3 tools have a solo_LTR hit
     - Intact clusters: all 3 tools have an intact_LTR_RT hit
  3. Family assignment: majority vote (2/3) on family name.
     If all three disagree: skip cluster (ambiguous family).
  4. For intact clusters: LTR boundary consensus — use the EDTA boundaries (most accurate
     from LTR_retriever) if available, otherwise average of available tools.
     Require LTR boundary agreement within MAX_BOUNDARY_DIFF bp (default 50).
  5. Assign final locus_id with species prefix + sequential number.

Output: catalog TSV with n_tools + per-tool vote columns added.

Usage:
    python scripts/semisyn/merge_annotations.py \\
        --edta       data/semisyn/annotations/rice/edta_loci.tsv \\
        --look4ltrs  data/semisyn/annotations/rice/look4ltrs_loci.tsv \\
        --soloLTRseeker data/semisyn/annotations/rice/soloLTRseeker_loci.tsv \\
        --species-prefix os \\
        --species rice \\
        --out data/semisyn/catalogs/rice_catalog.tsv
"""

import argparse
import csv
import sys
from pathlib import Path


MIN_OVERLAP = 0.50       # reciprocal overlap fraction required
MAX_BOUNDARY_DIFF = 50   # bp — maximum allowed LTR boundary deviation for intact elements

CATALOG_FIELDS = [
    "locus_id", "type", "family", "superfamily",
    "chrom", "start", "end", "strand",
    "ltr_5_start", "ltr_5_end", "ltr_3_start", "ltr_3_end",
    "int_start", "int_end",
    "n_tools", "edta", "look4ltrs", "soloLTRseeker",
]

TOOL_NAMES = ["edta", "look4ltrs", "soloLTRseeker"]


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def load_loci(path: Path, tool: str) -> list[dict]:
    rows = []
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            row["_tool"] = tool
            row["start"] = int(row["start"])
            row["end"]   = int(row["end"])
            for field in ("ltr_5_start", "ltr_5_end", "ltr_3_start", "ltr_3_end",
                          "int_start", "int_end"):
                val = row.get(field, ".")
                row[field] = int(val) if val not in (".", "", "NA") else None
            rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Overlap clustering
# ---------------------------------------------------------------------------

def reciprocal_overlap(a: dict, b: dict) -> float:
    """Reciprocal overlap between two loci (same chromosome assumed)."""
    ov = min(a["end"], b["end"]) - max(a["start"], b["start"]) + 1
    if ov <= 0:
        return 0.0
    return ov / min(a["end"] - a["start"] + 1, b["end"] - b["start"] + 1)


def cluster_loci(all_loci: list[dict]) -> list[list[dict]]:
    """
    Single-linkage clustering by reciprocal overlap >= MIN_OVERLAP.
    Only clusters loci of the same type (solo_LTR / intact_LTR_RT) on the same chrom.
    Returns list of clusters, each cluster is a list of loci.
    """
    # Group by (chrom, type) for efficiency
    from collections import defaultdict
    buckets: dict[tuple, list[dict]] = defaultdict(list)
    for locus in all_loci:
        buckets[(locus["chrom"], locus["type"])].append(locus)

    clusters = []
    for (chrom, ltype), loci in buckets.items():
        # Sort by start position for efficiency
        loci.sort(key=lambda x: x["start"])

        # Simple O(n^2) single-linkage clustering (adequate for gene-scale n)
        assigned = [-1] * len(loci)
        n_clusters = 0

        for i in range(len(loci)):
            if assigned[i] >= 0:
                continue
            assigned[i] = n_clusters
            n_clusters += 1
            # Grow cluster
            for j in range(i + 1, len(loci)):
                # Early exit: if loci[j].start is far past any current cluster member
                if loci[j]["start"] > loci[i]["end"] * 2:
                    break
                if assigned[j] >= 0:
                    # Already assigned — check if we can merge clusters
                    if reciprocal_overlap(loci[i], loci[j]) >= MIN_OVERLAP:
                        # Merge cluster j's cluster into cluster i's cluster
                        old_cid = assigned[j]
                        new_cid = assigned[i]
                        for k in range(len(loci)):
                            if assigned[k] == old_cid:
                                assigned[k] = new_cid
                else:
                    if reciprocal_overlap(loci[i], loci[j]) >= MIN_OVERLAP:
                        assigned[j] = assigned[i]

        # Collect clusters
        cluster_map: dict[int, list[dict]] = {}
        for i, cid in enumerate(assigned):
            cluster_map.setdefault(cid, []).append(loci[i])
        clusters.extend(cluster_map.values())

    return clusters


# ---------------------------------------------------------------------------
# Consensus derivation
# ---------------------------------------------------------------------------

def majority_vote(values: list[str]) -> str | None:
    """Return value with >1 vote, or None if all disagree."""
    from collections import Counter
    counts = Counter(values)
    top, top_count = counts.most_common(1)[0]
    return top if top_count > 1 else None


def consensus_coords(cluster: list[dict]) -> tuple[int, int, str]:
    """Return (start, end, strand) as median of cluster members."""
    starts = sorted(l["start"] for l in cluster)
    ends   = sorted(l["end"] for l in cluster)
    n = len(cluster)
    start = starts[n // 2]
    end   = ends[n // 2]
    strand = majority_vote([l["strand"] for l in cluster]) or "+"
    return start, end, strand


def consensus_boundaries(cluster: list[dict], prefer_edta: bool = True) -> dict | None:
    """
    Derive LTR boundary consensus for intact elements.

    Returns None if boundary agreement is too poor (> MAX_BOUNDARY_DIFF bp deviation).
    Prefers EDTA boundaries (most accurate, from LTR_retriever).
    """
    intact = [l for l in cluster if l["ltr_5_start"] is not None]
    if not intact:
        return None

    # Prefer EDTA if available
    if prefer_edta:
        edta_intact = [l for l in intact if l["_tool"] == "edta"]
        if edta_intact:
            ref = edta_intact[0]
        else:
            ref = intact[0]
    else:
        ref = intact[0]

    # Check boundary agreement across tools
    for l in intact:
        for field in ("ltr_5_start", "ltr_5_end", "ltr_3_start", "ltr_3_end"):
            if ref[field] is None or l[field] is None:
                continue
            if abs(ref[field] - l[field]) > MAX_BOUNDARY_DIFF:
                return None  # Too divergent — skip this cluster

    return {
        "ltr_5_start": ref["ltr_5_start"],
        "ltr_5_end":   ref["ltr_5_end"],
        "ltr_3_start": ref["ltr_3_start"],
        "ltr_3_end":   ref["ltr_3_end"],
        "int_start":   ref["int_start"],
        "int_end":     ref["int_end"],
    }


# ---------------------------------------------------------------------------
# Main merge
# ---------------------------------------------------------------------------

def merge(
    edta_path: Path | None,
    look4ltrs_path: Path | None,
    soloLTRseeker_path: Path | None,
    species_prefix: str,
    species: str,
    strict: bool = False,
) -> tuple[list[dict], dict]:
    """
    Returns (catalog_rows, stats_dict).
    """
    all_loci = []
    tools_present = []

    for path, tool in [
        (edta_path, "edta"),
        (look4ltrs_path, "look4ltrs"),
        (soloLTRseeker_path, "soloLTRseeker"),
    ]:
        if path is not None and path.exists():
            loci = load_loci(path, tool)
            all_loci.extend(loci)
            tools_present.append(tool)
            print(f"  Loaded {len(loci)} loci from {tool} ({path.name})")
        else:
            print(f"  Warning: {tool} loci not provided or file missing", file=sys.stderr)

    n_tools_available = len(tools_present)
    if n_tools_available < 3 and strict:
        print(f"ERROR: --strict mode requires all 3 tools; only {tools_present} provided.",
              file=sys.stderr)
        sys.exit(1)
    if n_tools_available < 2:
        print("ERROR: need at least 2 tool outputs to merge.", file=sys.stderr)
        sys.exit(1)

    print(f"\nClustering {len(all_loci)} total loci ...", flush=True)
    clusters = cluster_loci(all_loci)
    print(f"  {len(clusters)} clusters formed")

    catalog_rows = []
    solo_counter: dict[str, int] = {}
    intact_counter: dict[str, int] = {}

    stats = {
        "clusters_total": len(clusters),
        "clusters_all3": 0,
        "clusters_2tools": 0,
        "clusters_dropped_family": 0,
        "clusters_dropped_boundary": 0,
        "solo_accepted": 0,
        "intact_accepted": 0,
    }

    for cluster in clusters:
        cluster_tools = set(l["_tool"] for l in cluster)
        n = len(cluster_tools)

        if n < n_tools_available and strict:
            continue  # strict: require all available tools
        if n < 2:
            continue  # always require at least 2 tools

        if n == 3:
            stats["clusters_all3"] += 1
        else:
            stats["clusters_2tools"] += 1

        ltype = cluster[0]["type"]  # all same type by construction

        # Family: majority vote
        families = [l["family"] for l in cluster]
        family = majority_vote(families)
        if family is None:
            stats["clusters_dropped_family"] += 1
            continue

        # Superfamily: from family-matched locus
        superfamily = next(
            (l["superfamily"] for l in cluster if l["family"] == family),
            "unknown",
        )

        chrom = cluster[0]["chrom"]
        start, end, strand = consensus_coords(cluster)

        # Build vote flags
        tool_votes = {t: 1 if t in cluster_tools else 0 for t in TOOL_NAMES}

        if ltype == "intact_LTR_RT":
            bounds = consensus_boundaries(cluster)
            if bounds is None:
                stats["clusters_dropped_boundary"] += 1
                continue
            intact_counter[family] = intact_counter.get(family, 0) + 1
            locus_id = f"{species_prefix}_{family}_i{intact_counter[family]:04d}"
            row = {
                "locus_id": locus_id,
                "type": "intact_LTR_RT",
                "family": family,
                "superfamily": superfamily,
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                **bounds,
                "n_tools": n,
                **tool_votes,
            }
            stats["intact_accepted"] += 1

        else:  # solo_LTR
            solo_counter[family] = solo_counter.get(family, 0) + 1
            locus_id = f"{species_prefix}_{family}_s{solo_counter[family]:04d}"
            row = {
                "locus_id": locus_id,
                "type": "solo_LTR",
                "family": family,
                "superfamily": superfamily,
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "ltr_5_start": ".", "ltr_5_end": ".",
                "ltr_3_start": ".", "ltr_3_end": ".",
                "int_start": ".", "int_end": ".",
                "n_tools": n,
                **tool_votes,
            }
            stats["solo_accepted"] += 1

        catalog_rows.append(row)

    return catalog_rows, stats


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    global MIN_OVERLAP, MAX_BOUNDARY_DIFF
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--edta", type=Path, default=None,
                        help="EDTA per-tool loci TSV (from parse_edta.py)")
    parser.add_argument("--look4ltrs", type=Path, default=None,
                        help="Look4LTRs per-tool loci TSV (from parse_look4ltrs.py)")
    parser.add_argument("--soloLTRseeker", type=Path, default=None,
                        help="soloLTRseeker per-tool loci TSV (from parse_soloLTRseeker.py)")
    parser.add_argument("--species-prefix", default="sp",
                        help="Short species prefix for locus IDs")
    parser.add_argument("--species", default="unknown",
                        help="Species name for output metadata")
    parser.add_argument("--strict", action="store_true",
                        help="Require ALL provided tools to agree (drop 2-tool clusters)")
    parser.add_argument("--min-overlap", type=float, default=MIN_OVERLAP,
                        help=f"Minimum reciprocal overlap fraction (default: {MIN_OVERLAP})")
    parser.add_argument("--max-boundary-diff", type=int, default=MAX_BOUNDARY_DIFF,
                        help=f"Max LTR boundary deviation bp (default: {MAX_BOUNDARY_DIFF})")
    parser.add_argument("--out", required=True, type=Path,
                        help="Output catalog TSV")
    args = parser.parse_args()

    MIN_OVERLAP       = args.min_overlap
    MAX_BOUNDARY_DIFF = args.max_boundary_diff

    if not any([args.edta, args.look4ltrs, args.soloLTRseeker]):
        parser.error("Provide at least two tool TSVs (--edta, --look4ltrs, --soloLTRseeker)")

    catalog, stats = merge(
        edta_path=args.edta,
        look4ltrs_path=args.look4ltrs,
        soloLTRseeker_path=args.soloLTRseeker,
        species_prefix=args.species_prefix,
        species=args.species,
        strict=args.strict,
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CATALOG_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(catalog)

    # Summary
    print(f"\n=== Merge summary ===")
    print(f"  Clusters total:           {stats['clusters_total']}")
    print(f"  All-3-tool clusters:      {stats['clusters_all3']}")
    print(f"  2-tool clusters:          {stats['clusters_2tools']}")
    print(f"  Dropped (family ambig):   {stats['clusters_dropped_family']}")
    print(f"  Dropped (boundary diff):  {stats['clusters_dropped_boundary']}")
    print(f"  Accepted solo LTRs:       {stats['solo_accepted']}")
    print(f"  Accepted intact LTR-RTs:  {stats['intact_accepted']}")
    print(f"\nOutput: {args.out} ({len(catalog)} loci)")

    families = sorted(set(r["family"] for r in catalog))
    if families:
        print(f"\n{'Family':<25} {'Solo':>6} {'Intact':>8}")
        for fam in families:
            n_s = sum(1 for r in catalog if r["family"] == fam and r["type"] == "solo_LTR")
            n_i = sum(1 for r in catalog if r["family"] == fam and r["type"] == "intact_LTR_RT")
            print(f"  {fam:<23} {n_s:>6} {n_i:>8}")


if __name__ == "__main__":
    main()
