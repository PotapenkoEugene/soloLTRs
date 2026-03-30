"""
extract_sequences.py — Extract per-locus sequences from a genome FASTA.

For each locus in the catalog:
  - Extracts the full locus with FLANK bp of flanking context on each side
    (for virtual genome composition and read simulation)
  - For intact LTR-RTs: also extracts the LTR and internal sub-sequences
    (for consensus library building)

Uses samtools faidx for extraction (genome must be faidx-indexed or will be indexed here).

Outputs:
  loci_with_flanking.fa   — per-locus sequence with flanking (for compose_genome.py)
  extracted_ltrs.fa       — LTR sequences from intact elements (for consensus lib)
  extracted_internals.fa  — Internal sequences from intact elements (for consensus lib)
  loci_metadata.json      — Per-locus metadata: original coords, extracted window, family, type

Usage:
    python scripts/semisyn/extract_sequences.py \\
        --catalog data/semisyn/catalogs/t2_mock_catalog.tsv \\
        --genome data/sim/t2/genome/t2_genome.fa \\
        --flank 700 \\
        --out-dir data/semisyn/extracted/t2_mock
"""

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path


FLANK_DEFAULT = 700  # read_len(150) + insert_size(500) + margin(50)


def index_genome(genome: Path) -> None:
    """Build .fai index if not present."""
    fai = Path(str(genome) + ".fai")
    if not fai.exists():
        print(f"Indexing {genome} ...", flush=True)
        result = subprocess.run(
            ["samtools", "faidx", str(genome)],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            print(f"ERROR indexing genome: {result.stderr}", file=sys.stderr)
            sys.exit(1)


def get_chrom_lengths(genome: Path) -> dict[str, int]:
    """Read chromosome lengths from .fai index."""
    fai = Path(str(genome) + ".fai")
    lengths = {}
    with open(fai) as f:
        for line in f:
            parts = line.strip().split("\t")
            lengths[parts[0]] = int(parts[1])
    return lengths


def fetch_sequence(genome: Path, chrom: str, start: int, end: int) -> str:
    """
    Fetch sequence from genome using samtools faidx.
    Coordinates are 1-based inclusive (samtools convention).
    Returns uppercase sequence string (no newlines).
    """
    region = f"{chrom}:{start}-{end}"
    result = subprocess.run(
        ["samtools", "faidx", str(genome), region],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"samtools faidx failed for {region}: {result.stderr}")
    lines = result.stdout.strip().split("\n")
    # Skip the FASTA header line
    seq = "".join(lines[1:]).upper()
    return seq


def load_catalog(path: Path) -> list[dict]:
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def write_fasta(records: list[tuple[str, str]], path: Path) -> None:
    """Write a list of (header, sequence) pairs as FASTA (60-char line wrap)."""
    with open(path, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", required=True, type=Path,
                        help="Consensus locus catalog TSV")
    parser.add_argument("--genome", required=True, type=Path,
                        help="Genome FASTA (will be faidx-indexed if needed)")
    parser.add_argument("--flank", type=int, default=FLANK_DEFAULT,
                        help=f"Flanking bp on each side (default: {FLANK_DEFAULT})")
    parser.add_argument("--out-dir", required=True, type=Path,
                        help="Output directory")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    index_genome(args.genome)
    chrom_lengths = get_chrom_lengths(args.genome)
    catalog = load_catalog(args.catalog)

    loci_flanking: list[tuple[str, str]] = []      # (header, seq) for loci_with_flanking.fa
    extracted_ltrs: list[tuple[str, str]] = []     # (header, seq) for extracted_ltrs.fa
    extracted_ints: list[tuple[str, str]] = []     # (header, seq) for extracted_internals.fa
    metadata: list[dict] = []

    skipped = 0
    for row in catalog:
        locus_id = row["locus_id"]
        ltype = row["type"]
        family = row["family"]
        superfamily = row["superfamily"]
        chrom = row["chrom"]
        start = int(row["start"])   # 1-based inclusive
        end = int(row["end"])       # 1-based inclusive
        strand = row["strand"]

        chrom_len = chrom_lengths.get(chrom)
        if chrom_len is None:
            print(f"Warning: chromosome '{chrom}' not in genome, skipping {locus_id}",
                  file=sys.stderr)
            skipped += 1
            continue

        # Flanked window (clamp to chromosome boundaries)
        win_start = max(1, start - args.flank)
        win_end = min(chrom_len, end + args.flank)

        try:
            win_seq = fetch_sequence(args.genome, chrom, win_start, win_end)
        except RuntimeError as e:
            print(f"Warning: {e}, skipping {locus_id}", file=sys.stderr)
            skipped += 1
            continue

        # Flanked sequence header encodes original coords for compose_genome.py
        flank_header = (
            f"{locus_id} type={ltype} family={family} superfamily={superfamily} "
            f"chrom={chrom} locus_start={start} locus_end={end} strand={strand} "
            f"win_start={win_start} win_end={win_end}"
        )
        loci_flanking.append((flank_header, win_seq))

        meta = {
            "locus_id": locus_id,
            "type": ltype,
            "family": family,
            "superfamily": superfamily,
            "chrom": chrom,
            "locus_start": start,
            "locus_end": end,
            "strand": strand,
            "win_start": win_start,
            "win_end": win_end,
            "win_len": win_end - win_start + 1,
        }

        # For intact elements: also extract LTR and internal sub-sequences
        if ltype == "intact_LTR_RT" and row["ltr_5_start"] != ".":
            ltr_5_start = int(row["ltr_5_start"])
            ltr_5_end = int(row["ltr_5_end"])
            ltr_3_start = int(row["ltr_3_start"])
            ltr_3_end = int(row["ltr_3_end"])
            int_start = int(row["int_start"])
            int_end = int(row["int_end"])

            try:
                ltr_seq = fetch_sequence(args.genome, chrom, ltr_5_start, ltr_5_end)
                int_seq = fetch_sequence(args.genome, chrom, int_start, int_end)
            except RuntimeError as e:
                print(f"Warning: sub-element extraction failed for {locus_id}: {e}",
                      file=sys.stderr)
            else:
                # Header format: family#LTR/superfamily (matching make_lib.py convention)
                ltr_header = f"{locus_id}_{family}_LTR family={family} superfamily={superfamily}"
                int_header = f"{locus_id}_{family}_INT family={family} superfamily={superfamily}"
                extracted_ltrs.append((ltr_header, ltr_seq))
                extracted_ints.append((int_header, int_seq))

                meta.update({
                    "ltr_5_start": ltr_5_start, "ltr_5_end": ltr_5_end,
                    "ltr_3_start": ltr_3_start, "ltr_3_end": ltr_3_end,
                    "int_start": int_start, "int_end": int_end,
                    "ltr_len": ltr_5_end - ltr_5_start + 1,
                    "int_len": int_end - int_start + 1,
                })

        metadata.append(meta)

    # Write outputs
    write_fasta(loci_flanking, args.out_dir / "loci_with_flanking.fa")
    write_fasta(extracted_ltrs, args.out_dir / "extracted_ltrs.fa")
    write_fasta(extracted_ints, args.out_dir / "extracted_internals.fa")
    (args.out_dir / "loci_metadata.json").write_text(
        json.dumps(metadata, indent=2)
    )

    n_solo = sum(1 for r in catalog if r["type"] == "solo_LTR")
    n_intact = sum(1 for r in catalog if r["type"] == "intact_LTR_RT")
    print(f"Catalog: {len(catalog)} loci ({n_solo} solo, {n_intact} intact)")
    print(f"Extracted: {len(loci_flanking)} flanked sequences, "
          f"{len(extracted_ltrs)} LTRs, {len(extracted_ints)} internals")
    if skipped:
        print(f"Skipped: {skipped} loci (chromosome not found or extraction error)")
    print(f"Output: {args.out_dir}/")


if __name__ == "__main__":
    main()
