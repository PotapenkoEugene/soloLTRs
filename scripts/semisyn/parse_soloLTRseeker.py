"""
parse_soloLTRseeker.py — Parse soloLTRseeker output into the standardized loci catalog TSV.

soloLTRseeker (Orozco-Arias et al., https://github.com/eib-lab/soloLTRseeker) uses a
BLAST-based approach to identify solo LTRs relative to a known intact LTR library.
It requires LTR_retriever (or similar) to have previously identified intact elements.

Expected soloLTRseeker input files:
  --solo-gff3      GFF3 output from soloLTRseeker (solo LTRs)
  --intact-gff3    GFF3 with intact LTR-RTs (typically from LTR_retriever, which is
                   also a dependency of EDTA; can reuse *.EDTA.intact.gff3)

GFF3 format (soloLTRseeker solo_LTR entries):
  Chr1  soloLTRseeker  solo_LTR  100  500  .  +  .  ID=...;Name=FamilyName;...

Note: soloLTRseeker only identifies solo LTRs. Intact LTR-RT entries come from a
separate tool (LTR_retriever / EDTA). Use --intact-gff3 to provide the intact catalog
in any of the following formats:
  - EDTA intact.gff3 (repeat_region > LTR_retrotransposon > long_terminal_repeat ×2)
  - LTR_retriever pass.list (tab-separated text) via --intact-format ltr_retriever
  - Generic GFF3 with LTR_retrotransposon + long_terminal_repeat features (default)

Usage:
    python scripts/semisyn/parse_soloLTRseeker.py \\
        --solo-gff3   data/semisyn/annotations/rice/soloLTRseeker_solos.gff3 \\
        --intact-gff3 data/semisyn/annotations/rice/rice.EDTA.intact.gff3 \\
        --intact-format edta \\
        --species-prefix os \\
        --out data/semisyn/annotations/rice/soloLTRseeker_loci.tsv
"""

import argparse
import csv
import sys
from pathlib import Path


CATALOG_FIELDS = [
    "locus_id", "type", "family", "superfamily",
    "chrom", "start", "end", "strand",
    "ltr_5_start", "ltr_5_end", "ltr_3_start", "ltr_3_end",
    "int_start", "int_end",
]

SOLO_FEATURES   = {"solo_LTR", "soloLTR", "solo_long_terminal_repeat"}
INTACT_FEATURES = {"LTR_retrotransposon", "LTR_element", "intact_LTR",
                   "Copia_LTR_retrotransposon", "Gypsy_LTR_retrotransposon"}
LTR_BOUNDARY_FEATURES = {"long_terminal_repeat", "LTR"}


# ---------------------------------------------------------------------------
# Shared GFF3 utilities
# ---------------------------------------------------------------------------

def parse_attrs(attr_str: str) -> dict[str, str]:
    attrs = {}
    for part in attr_str.strip().split(";"):
        part = part.strip()
        if "=" in part:
            k, _, v = part.partition("=")
            attrs[k.strip()] = v.strip()
    return attrs


def parse_name(attrs: dict) -> tuple[str, str]:
    """
    Extract (family, superfamily) from GFF3 attributes.

    soloLTRseeker output uses:
      Name=soloLTR_01  (sequential number — not the family name)
      Classification=LTR/Gypsy_LTR_retrotransposon

    For soloLTRseeker entries, family is derived from Classification since Name
    is just a counter. Family name from Classification: strip the '_LTR_retrotransposon'
    suffix and use the main type (Gypsy, Copia, etc.) as a stand-in for the family.
    The merge step's majority vote will use EDTA/Look4LTRs names for final assignment.
    """
    name = attrs.get("Name", "")
    classification = attrs.get("Classification", attrs.get("Class", ""))

    # soloLTRseeker: Name=soloLTR_NN is not a real family name
    is_counter = name.lower().startswith("sololtrs") or name.lower().startswith("sololtр")
    name_is_sequential = name.lower().startswith("sololtrs_") or (
        name.lower().startswith("sololtrs") and name.split("_")[-1].isdigit()
    )
    # Detect sequential Name (soloLTR_01, soloLTR_02, ...) — not a family name
    if "_" in name and name.split("_")[-1].isdigit():
        name_is_sequential = True

    if "#" in name and not name_is_sequential:
        family, _, cls = name.partition("#")
        superfamily = cls.split("/")[-1] if "/" in cls else cls
    elif name and not name_is_sequential:
        family = name
        superfamily = classification.split("/")[-1] if "/" in classification else "unknown"
    elif classification:
        # Derive from Classification=LTR/Gypsy_LTR_retrotransposon -> superfamily=Gypsy
        cls_part = classification.split("/")[-1] if "/" in classification else classification
        # Strip _LTR_retrotransposon, _retrotransposon, etc.
        for sfx in ("_LTR_retrotransposon", "_retrotransposon", "_LTR"):
            if cls_part.endswith(sfx):
                cls_part = cls_part[: -len(sfx)]
                break
        superfamily = cls_part
        family = cls_part  # best available name; merge step will replace with EDTA/Look4LTRs name
    else:
        family = attrs.get("ID", "unknown")
        superfamily = "unknown"

    for suffix in ("-I", "_I", "-LTR", "_LTR", "#LTR", "-int", "_int"):
        if family.upper().endswith(suffix.upper()):
            family = family[: -len(suffix)]
            break

    return family, superfamily


def normalise_strand(s: str) -> str:
    return s if s in ("+", "-") else "+"


def iter_gff3(path: Path):
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, phase, attr_str = parts[:9]
            attrs = parse_attrs(attr_str)
            yield chrom, source, feature, int(start), int(end), strand, attrs


# ---------------------------------------------------------------------------
# Solo LTR parser
# ---------------------------------------------------------------------------

def parse_solo_gff3(path: Path, species_prefix: str) -> list[dict]:
    """Parse soloLTRseeker solo LTR GFF3."""
    rows = []
    solo_counter: dict[str, int] = {}

    for chrom, source, feature, start, end, strand, attrs in iter_gff3(path):
        # Accept solo_LTR and also bare LTR entries (some versions use different names)
        if feature not in SOLO_FEATURES and feature != "LTR":
            # If the source is soloLTRseeker, include all features as solo candidates
            if source.lower() not in ("soloLTRseeker", "sololtrs", "solo_ltr_seeker"):
                continue

        family, superfamily = parse_name(attrs)
        solo_counter[family] = solo_counter.get(family, 0) + 1
        locus_id = f"{species_prefix}_{family}_s{solo_counter[family]:04d}"
        rows.append({
            "locus_id": locus_id,
            "type": "solo_LTR",
            "family": family,
            "superfamily": superfamily,
            "chrom": chrom,
            "start": start,
            "end": end,
            "strand": normalise_strand(strand),
            "ltr_5_start": ".",
            "ltr_5_end":   ".",
            "ltr_3_start": ".",
            "ltr_3_end":   ".",
            "int_start":   ".",
            "int_end":     ".",
        })

    return rows


# ---------------------------------------------------------------------------
# Intact parser (generic GFF3 — same logic as EDTA for LTR_retrotransposon hierarchy)
# ---------------------------------------------------------------------------

def parse_intact_gff3(path: Path, species_prefix: str) -> list[dict]:
    """
    Parse intact LTR-RT GFF3 (accepts EDTA or LTR_retriever-style GFF3).
    Delegates to parent-child hierarchy parsing (same as parse_edta.py).
    """
    by_id: dict[str, dict] = {}
    children: dict[str, list[str]] = {}

    for chrom, source, feature, start, end, strand, attrs in iter_gff3(path):
        eid = attrs.get("ID", "")
        parent = attrs.get("Parent", "")
        if eid:
            by_id[eid] = {
                "chrom": chrom, "feature": feature,
                "start": start, "end": end,
                "strand": strand, "attrs": attrs,
            }
        if parent:
            children.setdefault(parent, []).append(eid)

    rows = []
    intact_counter: dict[str, int] = {}

    for eid, entry in by_id.items():
        if entry["feature"] not in INTACT_FEATURES:
            continue

        attrs = entry["attrs"]
        family, superfamily = parse_name(attrs)

        # In EDTA v2.0.0, LTRs are siblings under the same repeat_region parent
        parent_id = entry["attrs"].get("Parent", "")
        ltr_kids = []
        for cid in children.get(parent_id, []):
            child = by_id.get(cid, {})
            if child.get("feature") in LTR_BOUNDARY_FEATURES:
                ltr_kids.append(child)
        # Fallback: direct children (older EDTA versions with nested hierarchy)
        if not ltr_kids:
            for cid in children.get(eid, []):
                child = by_id.get(cid, {})
                if child.get("feature") in LTR_BOUNDARY_FEATURES:
                    ltr_kids.append(child)

        if len(ltr_kids) < 2:
            continue

        ltr_kids.sort(key=lambda x: x["start"])
        ltr5, ltr3 = ltr_kids[0], ltr_kids[-1]
        int_start = ltr5["end"] + 1
        int_end   = ltr3["start"] - 1
        if int_end <= int_start:
            continue

        intact_counter[family] = intact_counter.get(family, 0) + 1
        locus_id = f"{species_prefix}_{family}_i{intact_counter[family]:04d}"

        rows.append({
            "locus_id": locus_id,
            "type": "intact_LTR_RT",
            "family": family,
            "superfamily": superfamily,
            "chrom": entry["chrom"],
            "start": entry["start"],
            "end": entry["end"],
            "strand": normalise_strand(entry["strand"]),
            "ltr_5_start": ltr5["start"],
            "ltr_5_end":   ltr5["end"],
            "ltr_3_start": ltr3["start"],
            "ltr_3_end":   ltr3["end"],
            "int_start":   int_start,
            "int_end":     int_end,
        })

    return rows


# ---------------------------------------------------------------------------
# LTR_retriever pass.list parser (alternative intact format)
# ---------------------------------------------------------------------------

def parse_ltr_retriever_passlist(path: Path, species_prefix: str) -> list[dict]:
    """
    Parse LTR_retriever *.pass.list file (tab-separated, no header).

    Column order (LTR_retriever output):
      LTRid  chrom  ltr_start  ltr_end  strand  ltr_5_start  ltr_5_end
      ltr_3_start  ltr_3_end  int_start  int_end  family  ... (extra cols ignored)

    Note: coordinates are 1-based in LTR_retriever output.
    """
    rows = []
    intact_counter: dict[str, int] = {}

    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 11:
                continue
            try:
                chrom      = parts[1]
                start      = int(parts[2])
                end        = int(parts[3])
                strand     = normalise_strand(parts[4])
                ltr_5_start = int(parts[5])
                ltr_5_end   = int(parts[6])
                ltr_3_start = int(parts[7])
                ltr_3_end   = int(parts[8])
                int_start   = int(parts[9])
                int_end     = int(parts[10])
                family_raw  = parts[11] if len(parts) > 11 else "unknown"
            except (ValueError, IndexError):
                continue

            # Parse family and superfamily from name like "RLG_Dasheng#LTR/Gypsy"
            if "#" in family_raw:
                family, _, cls = family_raw.partition("#")
                superfamily = cls.split("/")[-1] if "/" in cls else cls
            else:
                family = family_raw
                superfamily = "unknown"

            intact_counter[family] = intact_counter.get(family, 0) + 1
            locus_id = f"{species_prefix}_{family}_i{intact_counter[family]:04d}"

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
                "ltr_5_end":   ltr_5_end,
                "ltr_3_start": ltr_3_start,
                "ltr_3_end":   ltr_3_end,
                "int_start":   int_start,
                "int_end":     int_end,
            })

    return rows


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solo-gff3", required=True, type=Path,
                        help="soloLTRseeker solo LTR GFF3 output")
    parser.add_argument("--intact-gff3", required=True, type=Path,
                        help="Intact LTR-RT GFF3 (EDTA or LTR_retriever format)")
    parser.add_argument("--intact-format", default="edta",
                        choices=["edta", "ltr_retriever"],
                        help="Format of the intact input (default: edta)")
    parser.add_argument("--species-prefix", default="sp",
                        help="Short species prefix for locus IDs")
    parser.add_argument("--out", required=True, type=Path,
                        help="Output loci TSV")
    args = parser.parse_args()

    print(f"Parsing solo LTRs from {args.solo_gff3} ...", flush=True)
    solo_rows = parse_solo_gff3(args.solo_gff3, args.species_prefix)
    print(f"  {len(solo_rows)} solo LTR loci")

    print(f"Parsing intact elements from {args.intact_gff3} ({args.intact_format}) ...",
          flush=True)
    if args.intact_format == "ltr_retriever":
        intact_rows = parse_ltr_retriever_passlist(args.intact_gff3, args.species_prefix)
    else:
        intact_rows = parse_intact_gff3(args.intact_gff3, args.species_prefix)
    print(f"  {len(intact_rows)} intact LTR-RT loci")

    all_rows = solo_rows + intact_rows

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CATALOG_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(all_rows)

    families = sorted(set(r["family"] for r in all_rows))
    print(f"\nWritten {len(all_rows)} loci to {args.out}")
    print(f"{'Family':<25} {'Solo':>6} {'Intact':>8}")
    for fam in families:
        n_s = sum(1 for r in all_rows if r["family"] == fam and r["type"] == "solo_LTR")
        n_i = sum(1 for r in all_rows if r["family"] == fam and r["type"] == "intact_LTR_RT")
        print(f"  {fam:<23} {n_s:>6} {n_i:>8}")


if __name__ == "__main__":
    main()
