"""
Tier 1: Minimal synthetic genome (~1 Mb, single chromosome).

Purpose: Baseline sensitivity — the simplest real-size genome. If a tool
         can't find solo LTRs here, it is fundamentally broken.

Genome specs:
  - 1 Mb background, 42% GC, 1 chromosome
  - 50 solo LTRs from RIRE2 (Ty3/Gypsy, 550 bp)
    - 10 at 0-1% divergence
    - 10 at 1-3% divergence
    - 10 at 3-5% divergence
    - 10 at 5-10% divergence
    - 10 at 10-15% divergence
  - 5 intact LTR-RTs from RIRE2
  - No other TEs

Outputs (data/sim/t1/):
  genome/t1_genome.fa           — gitignored
  ground_truth/t1_ground_truth.gff3
  ground_truth/t1_ground_truth.bed
  ground_truth/t1_ground_truth.json

Read simulation (ART, separate Makefile target):
  art_illumina -ss HS25 -i t1_genome.fa -l 150 -f 30 -m 500 -s 50 -p -o t1_reads
"""

import argparse
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[2]))

from scripts.sim.utils import (
    GenomeBuilder,
    generate_background,
    mutate_sequence,
    make_tsd,
    write_fasta,
    write_gff3,
    write_bed,
    write_manifest,
)


OUTDIR_DEFAULT = "data/sim/t1"
SEED = 100
GENOME_SIZE = 1_000_000


# Divergence groups: (div_min, div_max, count)
DIVERGENCE_GROUPS = [
    (0.00, 0.01, 10),
    (0.01, 0.03, 10),
    (0.03, 0.05, 10),
    (0.05, 0.10, 10),
    (0.10, 0.15, 10),
]

N_INTACT = 5


def _load_ltr(lib_path: Path, name: str) -> str:
    seqs: dict[str, list[str]] = {}
    current = None
    with open(lib_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                current = line[1:].split("#")[0].split()[0]
                seqs[current] = []
            elif current:
                seqs[current].append(line)
    if name not in seqs:
        raise KeyError(f"'{name}' not found in {lib_path}")
    return "".join(seqs[name])


def _load_internal(lib_path: Path, name: str) -> str:
    seqs: dict[str, list[str]] = {}
    current = None
    with open(lib_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                key = line[1:].split("#")[0].split("_INT")[0].split()[0]
                current = key
                seqs[key] = []
            elif current:
                seqs[current].append(line)
    if name not in seqs:
        raise KeyError(f"'{name}' not found in {lib_path}")
    return "".join(seqs[name])


def make_t1(lib: Path, outdir: Path) -> None:
    rng = random.Random(SEED)

    ltr = _load_ltr(lib / "ltr_only.fa", "RIRE2_LTR")
    internal = _load_internal(lib / "internal_only.fa", "RIRE2")

    bg = generate_background(GENOME_SIZE, gc_content=0.42, seed=SEED)

    # Build list of (divergence, position) for all elements
    # Total elements: 50 solo + 5 intact = 55
    # Minimum spacing: ensure no overlap (LTR ~550 bp, intact ~7 kb)
    # Space available: 1 Mb / 55 ≈ 18 kb per slot — plenty of room
    n_elements = sum(c for _, _, c in DIVERGENCE_GROUPS) + N_INTACT
    slot_size = GENOME_SIZE // (n_elements + 2)
    positions = sorted(rng.sample(range(slot_size, GENOME_SIZE - slot_size, slot_size // 2), n_elements))

    builder = GenomeBuilder(bg, chrom="chr1", tier=1)
    pos_idx = 0

    # Insert solo LTRs
    elem_seed = SEED + 1000
    for div_min, div_max, count in DIVERGENCE_GROUPS:
        for i in range(count):
            div = rng.uniform(div_min, div_max)
            ltr_mut = mutate_sequence(ltr, div, seed=elem_seed)
            tsd = make_tsd(5, seed=elem_seed)
            strand = rng.choice(["+", "-"])
            builder.insert_solo_ltr(
                position=positions[pos_idx],
                ltr_seq=ltr_mut,
                family="RIRE2",
                superfamily="Ty3_Gypsy",
                strand=strand,
                divergence=div,
                tsd=tsd,
                tsd_status="perfect",
                tsd_3=tsd,
                tier=1,
            )
            pos_idx += 1
            elem_seed += 1

    # Insert intact LTR-RTs (very low divergence — 1-2% LTR divergence for age)
    for i in range(N_INTACT):
        tsd = make_tsd(5, seed=elem_seed)
        div = rng.uniform(0.01, 0.03)
        strand = rng.choice(["+", "-"])
        builder.insert_intact_ltr_rt(
            position=positions[pos_idx],
            ltr_seq=ltr,
            internal_seq=internal,
            family="RIRE2",
            superfamily="Ty3_Gypsy",
            strand=strand,
            divergence=div,
            tsd=tsd,
            tsd_status="perfect",
        )
        pos_idx += 1
        elem_seed += 1

    genome_seq = builder.get_sequence()
    elements = builder.get_elements()

    # Write genome
    genome_dir = outdir / "genome"
    genome_dir.mkdir(parents=True, exist_ok=True)
    write_fasta([("chr1 t1_genome", genome_seq)], genome_dir / "t1_genome.fa")

    # Write ground truth
    gt_dir = outdir / "ground_truth"
    gt_dir.mkdir(parents=True, exist_ok=True)
    write_gff3(elements, gt_dir / "t1_ground_truth.gff3")
    write_bed(elements, gt_dir / "t1_ground_truth.bed")
    write_manifest(elements, gt_dir / "t1_ground_truth.json")

    solos = [e for e in elements if e.elem_type == "solo_LTR"]
    intact = [e for e in elements if e.elem_type == "intact_LTR_RT"]
    print(f"Tier 1 genome: {len(genome_seq):,} bp")
    print(f"  Solo LTRs:      {len(solos)} (expected 50)")
    print(f"  Intact LTR-RTs: {len(intact)} (expected 5)")
    print(f"  Genome FASTA:   {genome_dir}/t1_genome.fa")
    print(f"  Ground truth:   {gt_dir}/")

    # Divergence distribution check
    from collections import Counter
    bin_counts = Counter(e.divergence_bin for e in solos)
    print("\n  Solo LTR counts by divergence bin:")
    for b in sorted(bin_counts):
        print(f"    {b}: {bin_counts[b]}")


def main():
    parser = argparse.ArgumentParser(description="Generate Tier 1 synthetic genome")
    parser.add_argument("--libdir", default="data/lib")
    parser.add_argument("--outdir", default=OUTDIR_DEFAULT)
    args = parser.parse_args()
    make_t1(Path(args.libdir), Path(args.outdir))


if __name__ == "__main__":
    main()
