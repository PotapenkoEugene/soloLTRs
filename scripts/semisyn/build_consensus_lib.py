"""
build_consensus_lib.py — Build a per-species consensus library for sololtrs prepare.

Strategy: use EDTA's pre-built consensus library directly (*.EDTA.TElib.fa), which
already contains high-quality LTR and internal consensuses derived by LTR_retriever.
This avoids redundant MSA and is appropriate when EDTA was run with a real genome.

The EDTA library uses sequence names in the format:
  FamilyName#LTR/Superfamily     (LTR consensus sequences)
  FamilyName#LTR/Superfamily-I   (internal region consensuses)  [some versions]
  FamilyName-INT#LTR/Superfamily (alternative internal naming)

This script:
  1. Reads *.EDTA.TElib.fa
  2. Splits entries into LTR sequences vs. internal sequences by name
  3. Filters: keep only families present in the loci catalog (--catalog) with
     >= MIN_INTACT families (i.e., families with real ground truth)
  4. Writes separate ltr_only.fa and internal_only.fa compatible with `sololtrs prepare`

The output is equivalent to what `sololtrs prepare --lib` expects:
  ltr_only.fa      — one consensus LTR sequence per family
  internal_only.fa — one consensus internal sequence per family

Usage:
    python scripts/semisyn/build_consensus_lib.py \\
        --edta-lib  data/semisyn/annotations/rice/rice.EDTA.TElib.fa \\
        --catalog   data/semisyn/catalogs/rice_catalog.tsv \\
        --out-dir   data/semisyn/lib/rice \\
        [--min-intact 1]
"""

import argparse
import csv
import sys
from pathlib import Path


MIN_INTACT_DEFAULT = 1   # keep families with at least this many intact elements in catalog


def load_fasta(path: Path) -> list[tuple[str, str]]:
    """Load FASTA as list of (header, sequence) tuples (preserving order)."""
    records = []
    current_header = None
    current_parts: list[str] = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if current_header is not None:
                    records.append((current_header, "".join(current_parts)))
                current_header = line[1:]
                current_parts = []
            else:
                current_parts.append(line.upper())
    if current_header is not None:
        records.append((current_header, "".join(current_parts)))
    return records


def write_fasta(records: list[tuple[str, str]], path: Path, wrap: int = 60) -> None:
    with open(path, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            for i in range(0, len(seq), wrap):
                f.write(seq[i:i+wrap] + "\n")


def parse_family_from_header(header: str) -> tuple[str, bool]:
    """
    Parse family name and whether it is an internal-region sequence.

    EDTA/RepBase naming conventions handled:
      RLC_Angela#LTR/Copia         → family=RLC_Angela, is_internal=False
      RLC_Angela-I#LTR/Copia       → family=RLC_Angela, is_internal=True
      RLC_Angela-INT#LTR/Copia     → family=RLC_Angela, is_internal=True
      RLC_Angela_INT#LTR/Copia     → family=RLC_Angela, is_internal=True
      RLG_Dasheng#LTR/Gypsy-I      → family=RLG_Dasheng, is_internal=True (suffix in class)
    """
    # Take first whitespace-delimited token as the ID
    name = header.split()[0]

    # Split on '#' to get the base name vs classification
    if "#" in name:
        base, cls = name.split("#", 1)
    else:
        base = name
        cls = ""

    # Check for internal-region suffixes in base name
    is_internal = False
    INTERNAL_SUFFIXES = ("-I", "_I", "-INT", "_INT", "-int", "_int")
    for suf in INTERNAL_SUFFIXES:
        if base.upper().endswith(suf.upper()):
            base = base[: -len(suf)]
            is_internal = True
            break

    # Also check if the classification explicitly marks it as internal
    if not is_internal and cls.endswith("-I"):
        is_internal = True

    return base, is_internal


def load_catalog_families(path: Path, min_intact: int) -> set[str]:
    """Return set of family names with >= min_intact intact elements in catalog."""
    from collections import defaultdict
    counts: dict[str, int] = defaultdict(int)
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["type"] == "intact_LTR_RT":
                counts[row["family"]] += 1
    return {fam for fam, n in counts.items() if n >= min_intact}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--edta-lib", required=True, type=Path,
                        help="EDTA consensus library FASTA (*.EDTA.TElib.fa)")
    parser.add_argument("--catalog", required=True, type=Path,
                        help="Loci catalog TSV (from merge_annotations.py)")
    parser.add_argument("--out-dir", required=True, type=Path,
                        help="Output directory for ltr_only.fa and internal_only.fa")
    parser.add_argument("--min-intact", type=int, default=MIN_INTACT_DEFAULT,
                        help=f"Min intact elements per family to include (default: {MIN_INTACT_DEFAULT})")
    parser.add_argument("--no-filter", action="store_true",
                        help="Include all LTR families from the library (skip catalog filter)")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading EDTA library from {args.edta_lib} ...", flush=True)
    all_records = load_fasta(args.edta_lib)
    print(f"  {len(all_records)} sequences in library")

    # Load catalog families (for filtering)
    if not args.no_filter:
        allowed_families = load_catalog_families(args.catalog, args.min_intact)
        print(f"  Catalog families with >= {args.min_intact} intact elements: "
              f"{len(allowed_families)}")
    else:
        allowed_families = None

    # Split into LTR vs internal, filter to LTR class only
    ltr_records: dict[str, tuple[str, str]] = {}     # family -> (header, seq)
    internal_records: dict[str, tuple[str, str]] = {}

    skipped_non_ltr = 0
    skipped_not_in_catalog = 0

    for header, seq in all_records:
        # Only include LTR retrotransposons
        if "#LTR/" not in header and "#LTR#" not in header:
            # Check if it could be an LTR-RT without the standard naming
            if not any(kw in header.upper() for kw in ("GYPSY", "COPIA", "RETROVIRUS")):
                skipped_non_ltr += 1
                continue

        family, is_internal = parse_family_from_header(header)

        # Filter by catalog.  When the catalog uses superfamily-level names
        # (e.g. Arabidopsis with EDTA internal IDs TE_XXXXXXX normalised to
        # Copia/Gypsy), fall back to matching the library entry's superfamily
        # against the catalog family set.
        if allowed_families is not None and family not in allowed_families:
            if "#" in header:
                cls_part = header.split("#", 1)[1].split()[0]  # e.g. "LTR/Copia"
                sfam = cls_part.split("/")[-1] if "/" in cls_part else cls_part
                for sfx in ("-I", "_I"):
                    if sfam.endswith(sfx):
                        sfam = sfam[: -len(sfx)]
                        break
                if sfam in allowed_families:
                    family = sfam  # remap to superfamily for grouping
                else:
                    skipped_not_in_catalog += 1
                    continue
            else:
                skipped_not_in_catalog += 1
                continue

        if is_internal:
            # Keep longest internal per family (some EDTA libs have duplicates)
            if family not in internal_records or len(seq) > len(internal_records[family][1]):
                internal_records[family] = (header, seq)
        else:
            # LTR consensus — keep longest per family
            if family not in ltr_records or len(seq) > len(ltr_records[family][1]):
                ltr_records[family] = (header, seq)

    # Families with both LTR and internal consensus
    complete_families = sorted(set(ltr_records) & set(internal_records))
    ltr_only_families = sorted(set(ltr_records) - set(internal_records))

    print(f"\nLibrary parsing:")
    print(f"  Skipped (not LTR class):    {skipped_non_ltr}")
    print(f"  Skipped (not in catalog):   {skipped_not_in_catalog}")
    print(f"  Complete (LTR + internal):  {len(complete_families)}")
    print(f"  LTR-only (no internal):     {len(ltr_only_families)}")

    if ltr_only_families:
        print(f"  Warning: {len(ltr_only_families)} families have no internal consensus — "
              "will be excluded from pipeline library unless --no-filter is used with "
              "a library that contains internal sequences separately.", file=sys.stderr)

    # Write output files
    ltr_out_records = [(h, s) for fam, (h, s) in sorted(ltr_records.items())
                       if fam in complete_families]
    int_out_records = [(h, s) for fam, (h, s) in sorted(internal_records.items())
                       if fam in complete_families]

    ltr_path = args.out_dir / "ltr_only.fa"
    int_path  = args.out_dir / "internal_only.fa"

    write_fasta(ltr_out_records, ltr_path)
    write_fasta(int_out_records, int_path)

    print(f"\nOutput:")
    print(f"  {ltr_path} — {len(ltr_out_records)} LTR consensus sequences")
    print(f"  {int_path} — {len(int_out_records)} internal consensus sequences")
    print(f"\nFamilies included ({len(complete_families)}):")
    for fam in complete_families:
        ltr_len = len(ltr_records[fam][1])
        int_len = len(internal_records[fam][1])
        print(f"  {fam:<30} LTR={ltr_len} bp  INT={int_len} bp")

    if not complete_families:
        print(
            "\nWARNING: no families written. Check that --catalog families match "
            "EDTA library family names (Name=FamilyName#LTR/... format).",
            file=sys.stderr,
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
