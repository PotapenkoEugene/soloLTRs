"""
parse_edta.py — Parse EDTA output GFF3 files into the standardized loci catalog TSV.

EDTA (Extensive de novo TE Annotator) produces two relevant output files when run
with `--anno 1`:

  *.EDTA.intact.gff3
      Structurally complete LTR-RTs with sub-element boundaries annotated by
      LTR_retriever. Hierarchy: repeat_region > {Copia,Gypsy,}LTR_retrotransposon >
      long_terminal_repeat (×2). Used for intact_LTR_RT loci.

      EDTA v2.0.0 uses 'Copia_LTR_retrotransposon' and 'Gypsy_LTR_retrotransposon'
      as the element feature type (instead of the generic 'LTR_retrotransposon').

  *.EDTA.TEanno.gff3
      Full genome homology annotation. EDTA v2.0.0 does NOT produce explicit
      'solo_LTR' features here. Solo LTRs are inferred as LTR-type entries
      (Copia_/Gypsy_/LTR_retrotransposon) that:
        (a) do not overlap any structurally intact element, and
        (b) are shorter than MAX_SOLO_LEN (solo LTRs < one full element length).
      Older EDTA versions that do emit 'solo_LTR' features are handled automatically.

GFF3 attribute parsing:
  Name=RLC_Angela#LTR/Copia  →  family=RLC_Angela, superfamily=Copia
  Name=RLG_Dasheng#LTR/Gypsy →  family=RLG_Dasheng, superfamily=Gypsy
  Name=TE_00000701            →  internal EDTA ID — normalised to superfamily
  Classification=LTR/Copia   →  superfamily=Copia (fallback)
  Strand is often '?' in EDTA output — normalised to '+' (unknown = unstranded).

Usage:
    python scripts/semisyn/parse_edta.py \\
        --intact-gff3 data/semisyn/annotations/arabidopsis/TAIR10.fa.mod.EDTA.intact.gff3 \\
        --anno-gff3   data/semisyn/annotations/arabidopsis/TAIR10.fa.mod.EDTA.TEanno.gff3 \\
        --species-prefix at \\
        --out data/semisyn/annotations/arabidopsis/edta_loci.tsv
"""

import argparse
import csv
import re
import sys
from pathlib import Path


CATALOG_FIELDS = [
    "locus_id", "type", "family", "superfamily",
    "chrom", "start", "end", "strand",
    "ltr_5_start", "ltr_5_end", "ltr_3_start", "ltr_3_end",
    "int_start", "int_end",
]

# EDTA v2.0.0 prefixes element feature type with superfamily name
INTACT_LTR_FEATURES = {
    "LTR_retrotransposon",
    "Copia_LTR_retrotransposon",
    "Gypsy_LTR_retrotransposon",
}

# Same set for TEanno LTR annotations (used to find solo LTR candidates)
ANNO_LTR_FEATURES = INTACT_LTR_FEATURES | {"solo_LTR"}

# EDTA internal IDs (TE_00000XXX) are not biological family names — use superfamily instead
_EDTA_INTERNAL_ID = re.compile(r"^TE_\d{5,}")


# ---------------------------------------------------------------------------
# GFF3 utilities
# ---------------------------------------------------------------------------

def parse_attrs(attr_str: str) -> dict[str, str]:
    """Parse GFF3 column 9 attribute string into a dict."""
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

    Priority:
      1. Name=FamilyName#LTR/Superfamily  (named family, e.g. RLC_Angela#LTR/Copia)
      2. Name=TE_XXXXXXX (EDTA internal ID) → use Classification for superfamily
      3. Name=FamilyName (superfamily from Classification)
      4. Alias or ID as fallback

    For species where EDTA assigns internal IDs (e.g. Arabidopsis), family is
    normalised to superfamily so that tools can be merged by a common name.
    """
    name = attrs.get("Name", "")
    classification = attrs.get("Classification", "")

    # EDTA internal IDs are not biological family names; use superfamily from Classification
    if _EDTA_INTERNAL_ID.match(name):
        sfam_raw = classification.split("/")[-1] if "/" in classification else "unknown"
        for sfx in ("_LTR_retrotransposon", "_retrotransposon", "_LTR"):
            if sfam_raw.endswith(sfx):
                sfam_raw = sfam_raw[: -len(sfx)]
                break
        sfam = sfam_raw if sfam_raw else "unknown"
        return sfam, sfam

    if "#" in name:
        family, _, cls = name.partition("#")
        superfamily = cls.split("/")[-1] if "/" in cls else cls
    elif name:
        family = name
        superfamily = classification.split("/")[-1] if "/" in classification else "unknown"
    else:
        family = attrs.get("Alias", attrs.get("ID", "unknown"))
        superfamily = classification.split("/")[-1] if "/" in classification else "unknown"

    # Strip common RepBase/EDTA suffixes that aren't part of the family name
    for suffix in ("-I", "_I", "-LTR", "_LTR", "#LTR", "-int", "_int"):
        if family.upper().endswith(suffix.upper()):
            family = family[: -len(suffix)]
            break

    return family, superfamily


def normalise_strand(s: str) -> str:
    return s if s in ("+", "-") else "+"


def iter_gff3(path: Path):
    """Yield (chrom, source, feature, start, end, strand, attrs_dict) for non-comment lines."""
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
# Intact LTR-RT parser
# ---------------------------------------------------------------------------

def parse_intact_gff3(path: Path, species_prefix: str) -> list[dict]:
    """
    Parse *.EDTA.intact.gff3 into intact_LTR_RT catalog rows.

    Expected hierarchy (EDTA v2.0.0):
      repeat_region  (parent, coords = element coords)
        Copia_LTR_retrotransposon  (or Gypsy_ or plain LTR_retrotransposon)
          long_terminal_repeat  ×2  (grandchildren, LTR coords)

    The internal region is derived: int_start = ltr_5_end+1, int_end = ltr_3_start-1.
    """
    # First pass: collect all entries keyed by their GFF ID
    by_id: dict[str, dict] = {}
    children: dict[str, list[str]] = {}  # parent_id -> [child_ids]

    for chrom, source, feature, start, end, strand, attrs in iter_gff3(path):
        eid = attrs.get("ID", "")
        parent = attrs.get("Parent", "")
        entry = {
            "chrom": chrom,
            "feature": feature,
            "start": start,
            "end": end,
            "strand": strand,
            "attrs": attrs,
        }
        if eid:
            by_id[eid] = entry
        if parent:
            children.setdefault(parent, []).append(eid)

    # Second pass: collect LTR_retrotransposon entries and their LTR children
    rows = []
    intact_counter: dict[str, int] = {}

    for eid, entry in by_id.items():
        if entry["feature"] not in INTACT_LTR_FEATURES:
            continue

        attrs = entry["attrs"]
        family, superfamily = parse_name(attrs)

        # Check classification is actually LTR
        cls = attrs.get("Classification", "")
        if cls and not cls.upper().startswith("LTR"):
            continue

        # In EDTA v2.0.0, LTRs are siblings of the element under the same repeat_region.
        # Collect long_terminal_repeat entries from the parent repeat_region's children.
        parent_id = attrs.get("Parent", "")
        ltr_kids = []
        for cid in children.get(parent_id, []):
            child = by_id.get(cid, {})
            if child.get("feature") == "long_terminal_repeat":
                ltr_kids.append(child)
        # Also check direct children (for EDTA versions with nested hierarchy)
        if not ltr_kids:
            for cid in children.get(eid, []):
                child = by_id.get(cid, {})
                if child.get("feature") == "long_terminal_repeat":
                    ltr_kids.append(child)

        if len(ltr_kids) < 2:
            # Skip elements without two LTR boundaries (incomplete annotation)
            continue

        # Sort by position: first = 5' LTR, second = 3' LTR
        ltr_kids.sort(key=lambda x: x["start"])
        ltr5 = ltr_kids[0]
        ltr3 = ltr_kids[-1]

        ltr_5_start = ltr5["start"]
        ltr_5_end   = ltr5["end"]
        ltr_3_start = ltr3["start"]
        ltr_3_end   = ltr3["end"]
        int_start   = ltr_5_end + 1
        int_end     = ltr_3_start - 1

        if int_end <= int_start:
            # Degenerate element: LTRs overlap (shouldn't happen in EDTA output)
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
            "ltr_5_start": ltr_5_start,
            "ltr_5_end":   ltr_5_end,
            "ltr_3_start": ltr_3_start,
            "ltr_3_end":   ltr_3_end,
            "int_start":   int_start,
            "int_end":     int_end,
        })

    return rows


# ---------------------------------------------------------------------------
# Solo LTR parser (from TEanno.gff3)
# ---------------------------------------------------------------------------

def _build_intact_index(intact_rows: list[dict]) -> dict[str, list[tuple[int, int]]]:
    """Build chrom → sorted list of (start, end) for intact element intervals."""
    index: dict[str, list[tuple[int, int]]] = {}
    for row in intact_rows:
        index.setdefault(row["chrom"], []).append((row["start"], row["end"]))
    for chrom in index:
        index[chrom].sort()
    return index


def _overlaps_intact(
    chrom: str,
    start: int,
    end: int,
    intact_index: dict[str, list[tuple[int, int]]],
    min_overlap_frac: float = 0.5,
) -> bool:
    """Return True if (chrom, start, end) overlaps an intact element by >= min_overlap_frac."""
    length = end - start + 1
    for istart, iend in intact_index.get(chrom, []):
        if iend < start:
            continue
        if istart > end:
            break
        overlap = min(end, iend) - max(start, istart) + 1
        if overlap / length >= min_overlap_frac:
            return True
    return False


def parse_solo_gff3(path: Path, species_prefix: str, intact_rows: list[dict]) -> list[dict]:
    """
    Parse solo LTR entries from *.EDTA.TEanno.gff3.

    EDTA v2.0.0 strategy (no 'solo_LTR' feature):
      Collect Copia_/Gypsy_/LTR_retrotransposon entries from TEanno that:
        1. Do NOT overlap (>= 50%) any structurally intact element
        2. Are shorter than MAX_SOLO_LEN (filters out full intact-length entries
           missed by LTR_retriever and very long fragment chains)
      Entries < MIN_SOLO_LEN are discarded as spurious short fragments.

    Older EDTA versions that produce explicit 'solo_LTR' features use those directly.
    """
    MIN_SOLO_LEN = 100   # bp — discard sub-LTR fragments
    MAX_SOLO_LEN = 3000  # bp — solo LTRs are typically 100–2500 bp; full elements > 4 kb

    intact_index = _build_intact_index(intact_rows)
    explicit_solos = []
    candidate_entries = []

    for chrom, source, feature, start, end, strand, attrs in iter_gff3(path):
        if feature == "solo_LTR":
            explicit_solos.append((chrom, source, feature, start, end, strand, attrs))
            continue

        if feature not in ANNO_LTR_FEATURES:
            continue

        cls = attrs.get("Classification", "")
        if not cls.upper().startswith("LTR"):
            continue

        length = end - start + 1
        if length < MIN_SOLO_LEN or length > MAX_SOLO_LEN:
            continue

        if _overlaps_intact(chrom, start, end, intact_index):
            continue

        candidate_entries.append((chrom, source, feature, start, end, strand, attrs))

    if explicit_solos:
        source_entries = explicit_solos
    else:
        print(
            f"Note: no solo_LTR features in {path.name} (EDTA v2.0.0 mode). "
            f"Using {len(candidate_entries)} TEanno LTR entries "
            f"(length {MIN_SOLO_LEN}–{MAX_SOLO_LEN} bp, not overlapping intacts).",
            file=sys.stderr,
        )
        source_entries = candidate_entries

    rows = []
    solo_counter: dict[str, int] = {}
    for chrom, source, feature, start, end, strand, attrs in source_entries:
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
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--intact-gff3", required=True, type=Path,
                        help="*.EDTA.intact.gff3")
    parser.add_argument("--anno-gff3", required=True, type=Path,
                        help="*.EDTA.TEanno.gff3 (for solo LTRs)")
    parser.add_argument("--species-prefix", default="sp",
                        help="Short species prefix for locus IDs (e.g. 'at' for Arabidopsis)")
    parser.add_argument("--out", required=True, type=Path,
                        help="Output loci TSV")
    args = parser.parse_args()

    print(f"Parsing intact elements from {args.intact_gff3.name} ...", flush=True)
    intact_rows = parse_intact_gff3(args.intact_gff3, args.species_prefix)
    print(f"  {len(intact_rows)} intact LTR-RT loci")

    print(f"Parsing solo LTRs from {args.anno_gff3.name} ...", flush=True)
    solo_rows = parse_solo_gff3(args.anno_gff3, args.species_prefix, intact_rows)
    print(f"  {len(solo_rows)} solo LTR loci")

    all_rows = solo_rows + intact_rows  # solos first (more numerous, better for inspection)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CATALOG_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(all_rows)

    # Summary
    families = sorted(set(r["family"] for r in all_rows))
    print(f"\nWritten {len(all_rows)} loci to {args.out}")
    print(f"{'Family':<25} {'Solo':>6} {'Intact':>8}")
    for fam in families:
        n_s = sum(1 for r in all_rows if r["family"] == fam and r["type"] == "solo_LTR")
        n_i = sum(1 for r in all_rows if r["family"] == fam and r["type"] == "intact_LTR_RT")
        print(f"  {fam:<23} {n_s:>6} {n_i:>8}")


if __name__ == "__main__":
    main()
