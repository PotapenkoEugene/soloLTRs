"""
parse_look4ltrs.py — Parse Look4LTRs output into the standardized loci catalog TSV.

Look4LTRs (Orozco-Arias et al. 2024, Nucleic Acids Research doi:10.1093/nar/gkae510)
is a C++ tool that identifies both intact LTR-RTs and solo LTRs directly from a
genome FASTA. It explicitly classifies solo LTRs without requiring RepeatMasker.

Expected Look4LTRs output directory structure:
  <results_dir>/
    <genome>_intact.gff3         — intact LTR-RT loci with sub-element boundaries
    <genome>_solos.gff3          — solo LTR loci
    <genome>_internals.gff3      — internal regions (optional; used when present)

GFF3 format for intact elements (Look4LTRs convention):
  Chr1  Look4LTRs  LTR_retrotransposon  100  5000  .  +  .  ID=...;Name=FamilyName;...
  Chr1  Look4LTRs  long_terminal_repeat  100  500  .  +  .  Parent=...
  Chr1  Look4LTRs  long_terminal_repeat  4500  5000  .  +  .  Parent=...

GFF3 format for solo LTRs:
  Chr1  Look4LTRs  solo_LTR  100  500  .  +  .  ID=...;Name=FamilyName;...

Family name conventions:
  Look4LTRs uses the family name from the reference library provided or infers it
  de novo. Name attribute typically: FamilyName or FamilyName#LTR/Superfamily.

Note: Look4LTRs output format may vary between versions. If your version produces
different feature types, adjust INTACT_FEATURE and SOLO_FEATURE below.

Usage:
    python scripts/semisyn/parse_look4ltrs.py \\
        --results-dir data/semisyn/annotations/rice/look4ltrs/ \\
        --species-prefix os \\
        --out data/semisyn/annotations/rice/look4ltrs_loci.tsv
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

# Feature type strings used by Look4LTRs — adjust if your version differs
INTACT_FEATURES = {"LTR_retrotransposon", "LTR_element", "intact_LTR"}
SOLO_FEATURES   = {"solo_LTR", "solo_long_terminal_repeat", "LTR_solo"}
LTR_BOUNDARY_FEATURES = {"long_terminal_repeat", "LTR_boundary", "LTR"}


# ---------------------------------------------------------------------------
# Shared GFF3 utilities (same as parse_edta.py)
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
    name = attrs.get("Name", "")
    classification = attrs.get("Classification", attrs.get("Class", ""))

    if "#" in name:
        family, _, cls = name.partition("#")
        superfamily = cls.split("/")[-1] if "/" in cls else cls
    elif name:
        family = name
        superfamily = classification.split("/")[-1] if "/" in classification else "unknown"
    else:
        family = attrs.get("ID", "unknown")
        superfamily = classification.split("/")[-1] if "/" in classification else "unknown"

    for suffix in ("-I", "_I", "-LTR", "_LTR", "#LTR", "-int", "_int"):
        if family.upper().endswith(suffix.upper()):
            family = family[: -len(suffix)]
            break

    return family, superfamily


def normalise_strand(s: str) -> str:
    return s if s in ("+", "-") else "+"


def find_gff3(directory: Path, pattern: str) -> Path | None:
    """Return the first file matching the glob pattern, or None."""
    matches = sorted(directory.glob(pattern))
    return matches[0] if matches else None


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
# Intact parser
# ---------------------------------------------------------------------------

def parse_intact_gff3(path: Path, species_prefix: str) -> list[dict]:
    """Parse Look4LTRs intact GFF3 using parent-child hierarchy."""
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

        # Collect LTR boundary children
        ltr_kids = []
        for cid in children.get(eid, []):
            child = by_id.get(cid, {})
            if child.get("feature") in LTR_BOUNDARY_FEATURES:
                ltr_kids.append(child)

        if len(ltr_kids) < 2:
            # Fallback: if LTRs not in hierarchy, check for separate LTR annotations
            # that overlap this element's coordinates (handled elsewhere)
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
# Solo parser
# ---------------------------------------------------------------------------

def parse_solo_gff3(path: Path, species_prefix: str) -> list[dict]:
    rows = []
    solo_counter: dict[str, int] = {}

    for chrom, source, feature, start, end, strand, attrs in iter_gff3(path):
        if feature not in SOLO_FEATURES:
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
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", required=True, type=Path,
                        help="Look4LTRs output directory")
    parser.add_argument("--intact-gff3", type=Path, default=None,
                        help="Intact GFF3 (auto-detected from --results-dir if omitted)")
    parser.add_argument("--solo-gff3", type=Path, default=None,
                        help="Solo LTR GFF3 (auto-detected from --results-dir if omitted)")
    parser.add_argument("--species-prefix", default="sp",
                        help="Short species prefix for locus IDs")
    parser.add_argument("--out", required=True, type=Path,
                        help="Output loci TSV")
    args = parser.parse_args()

    # Auto-detect GFF3 files if not specified
    intact_gff3 = args.intact_gff3
    solo_gff3   = args.solo_gff3

    if intact_gff3 is None:
        intact_gff3 = find_gff3(args.results_dir, "*intact*.gff3")
        if intact_gff3 is None:
            intact_gff3 = find_gff3(args.results_dir, "*LTR*.gff3")
        if intact_gff3 is None:
            print("ERROR: could not auto-detect intact GFF3; use --intact-gff3",
                  file=sys.stderr)
            sys.exit(1)
        print(f"Auto-detected intact GFF3: {intact_gff3}")

    if solo_gff3 is None:
        solo_gff3 = find_gff3(args.results_dir, "*solo*.gff3")
        if solo_gff3 is None:
            print("ERROR: could not auto-detect solo LTR GFF3; use --solo-gff3",
                  file=sys.stderr)
            sys.exit(1)
        print(f"Auto-detected solo GFF3: {solo_gff3}")

    print(f"Parsing intact elements from {intact_gff3} ...", flush=True)
    intact_rows = parse_intact_gff3(intact_gff3, args.species_prefix)
    print(f"  {len(intact_rows)} intact LTR-RT loci")

    print(f"Parsing solo LTRs from {solo_gff3} ...", flush=True)
    solo_rows = parse_solo_gff3(solo_gff3, args.species_prefix)
    print(f"  {len(solo_rows)} solo LTR loci")

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
