"""
Tier 0: Unit test sequences (~10 KB total, seconds to run).

Generates 4 subcategories, each as a separate FASTA + GFF3 + BED:
  t0_clean/        — 5 perfect solo LTRs, 0% divergence, perfect TSDs
  t0_diverged/     — 5 solo LTRs, 10-20% divergence, perfect TSDs
  t0_intact_nearby/ — 3 solo LTRs + 2 flanking intact LTR-RTs (discrimination test)
  t0_edge_cases/   — <200 bp LTR (CRM1), degraded TSD, partial remnant

All categories use RIRE2 (Ty3/Gypsy, 550 bp LTR) as the primary family,
except edge_cases which also uses CRM1 (350 bp) to test shorter LTRs.

Background: 10 kb random sequence (42% GC) for each test.
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
    degrade_tsd,
    write_fasta,
    write_gff3,
    write_bed,
    write_manifest,
)


OUTDIR_DEFAULT = "data/sim/t0"
SEED_BASE = 42


def _load_ltr(lib_path: Path, family: str) -> str:
    """Extract LTR sequence for a given family from ltr_only.fa."""
    seqs = {}
    current = None
    with open(lib_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                current = line[1:].split("#")[0].split()[0]
                seqs[current] = []
            elif current:
                seqs[current].append(line)
    if family not in seqs:
        raise KeyError(f"Family '{family}' not found in {lib_path}")
    return "".join(seqs[family])


def _load_internal(lib_path: Path, family: str) -> str:
    seqs = {}
    current = None
    with open(lib_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                current = line[1:].split("#")[0].split("_INT")[0].split()[0]
                seqs[current] = []
            elif current:
                seqs[current].append(line)
    if family not in seqs:
        raise KeyError(f"Family '{family}' not found in {lib_path}")
    return "".join(seqs[family])


def make_t0_clean(ltr_rire2: str, outdir: Path, tier: int = 0) -> None:
    """5 perfect solo LTRs at 0% divergence with perfect TSDs in 10 kb background."""
    outdir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(SEED_BASE + 1)
    bg_len = 10_000
    bg = generate_background(bg_len, gc_content=0.42, seed=SEED_BASE + 1)
    builder = GenomeBuilder(bg, chrom="chr1", tier=tier)

    # Place 5 solo LTRs at evenly spaced positions (avoid 0 and end)
    positions = [1200, 2600, 4200, 6000, 7800]
    for i, pos in enumerate(positions):
        tsd = make_tsd(5, seed=SEED_BASE + 100 + i)
        builder.insert_solo_ltr(
            position=pos,
            ltr_seq=ltr_rire2,
            family="RIRE2",
            superfamily="Ty3_Gypsy",
            strand=rng.choice(["+", "-"]),
            divergence=0.0,
            tsd=tsd,
            tsd_status="perfect",
            tsd_3=tsd,
            tier=tier,
        )

    seq = builder.get_sequence()
    elements = builder.get_elements()

    write_fasta([("chr1 t0_clean", seq)], outdir / "t0_clean.fa")
    write_gff3(elements, outdir / "t0_clean.gff3")
    write_bed(elements, outdir / "t0_clean.bed")
    write_manifest(elements, outdir / "t0_clean.json")

    print(f"  t0_clean: {len(elements)} solo LTRs in {len(seq):,} bp genome")


def make_t0_diverged(ltr_rire2: str, outdir: Path, tier: int = 0) -> None:
    """5 solo LTRs at 10-20% divergence with perfect TSDs in 10 kb background."""
    outdir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(SEED_BASE + 2)
    bg = generate_background(10_000, gc_content=0.42, seed=SEED_BASE + 2)
    builder = GenomeBuilder(bg, chrom="chr1", tier=tier)

    divergences = [0.10, 0.12, 0.15, 0.17, 0.20]
    positions = [1000, 2800, 4500, 6200, 8000]
    for i, (pos, div) in enumerate(zip(positions, divergences)):
        tsd = make_tsd(5, seed=SEED_BASE + 200 + i)
        ltr_mut = mutate_sequence(ltr_rire2, div, seed=SEED_BASE + 200 + i)
        builder.insert_solo_ltr(
            position=pos,
            ltr_seq=ltr_mut,
            family="RIRE2",
            superfamily="Ty3_Gypsy",
            strand=rng.choice(["+", "-"]),
            divergence=div,
            tsd=tsd,
            tsd_status="perfect",
            tsd_3=tsd,
            tier=tier,
        )

    seq = builder.get_sequence()
    elements = builder.get_elements()

    write_fasta([("chr1 t0_diverged", seq)], outdir / "t0_diverged.fa")
    write_gff3(elements, outdir / "t0_diverged.gff3")
    write_bed(elements, outdir / "t0_diverged.bed")
    write_manifest(elements, outdir / "t0_diverged.json")

    print(f"  t0_diverged: {len(elements)} solo LTRs in {len(seq):,} bp genome")


def make_t0_intact_nearby(
    ltr_rire2: str,
    internal_rire2: str,
    outdir: Path,
    tier: int = 0,
) -> None:
    """3 solo LTRs + 2 intact LTR-RTs — discrimination test."""
    outdir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(SEED_BASE + 3)
    bg = generate_background(30_000, gc_content=0.42, seed=SEED_BASE + 3)
    builder = GenomeBuilder(bg, chrom="chr1", tier=tier)

    # 2 intact LTR-RTs first (at low positions)
    for i, pos in enumerate([500, 15000]):
        tsd = make_tsd(5, seed=SEED_BASE + 300 + i)
        builder.insert_intact_ltr_rt(
            position=pos,
            ltr_seq=ltr_rire2,
            internal_seq=internal_rire2,
            family="RIRE2",
            superfamily="Ty3_Gypsy",
            strand=rng.choice(["+", "-"]),
            divergence=0.02,
            tsd=tsd,
            tsd_status="perfect",
        )

    # 3 solo LTRs interspersed
    for i, pos in enumerate([8000, 12000, 22000]):
        tsd = make_tsd(5, seed=SEED_BASE + 310 + i)
        builder.insert_solo_ltr(
            position=pos,
            ltr_seq=ltr_rire2,
            family="RIRE2",
            superfamily="Ty3_Gypsy",
            strand=rng.choice(["+", "-"]),
            divergence=0.03,
            tsd=tsd,
            tsd_status="perfect",
            tsd_3=tsd,
            tier=tier,
        )

    seq = builder.get_sequence()
    elements = builder.get_elements()

    write_fasta([("chr1 t0_intact_nearby", seq)], outdir / "t0_intact_nearby.fa")
    write_gff3(elements, outdir / "t0_intact_nearby.gff3")
    write_bed(elements, outdir / "t0_intact_nearby.bed")
    write_manifest(elements, outdir / "t0_intact_nearby.json")

    solos = [e for e in elements if e.elem_type == "solo_LTR"]
    intact = [e for e in elements if e.elem_type == "intact_LTR_RT"]
    print(f"  t0_intact_nearby: {len(solos)} solo LTRs + {len(intact)} intact elements in {len(seq):,} bp genome")


def make_t0_edge_cases(ltr_rire2: str, ltr_crm1: str, outdir: Path, tier: int = 0) -> None:
    """Edge cases: short LTR, degraded TSD, partial internal remnant."""
    outdir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(SEED_BASE + 4)
    bg = generate_background(10_000, gc_content=0.42, seed=SEED_BASE + 4)
    builder = GenomeBuilder(bg, chrom="chr1", tier=tier)

    # Case 1: Short LTR (CRM1, ~350 bp)
    tsd1 = make_tsd(5, seed=SEED_BASE + 401)
    builder.insert_solo_ltr(
        position=800,
        ltr_seq=ltr_crm1,
        family="CRM1",
        superfamily="Ty3_Gypsy",
        strand="+",
        divergence=0.05,
        tsd=tsd1,
        tsd_status="perfect",
        tsd_3=tsd1,
        tier=tier,
    )

    # Case 2: 1-mismatch TSD
    raw_tsd = make_tsd(5, seed=SEED_BASE + 402)
    tsd_5, tsd_3 = degrade_tsd(raw_tsd, mode="mismatch_1", seed=SEED_BASE + 402)
    builder.insert_solo_ltr(
        position=2500,
        ltr_seq=ltr_rire2,
        family="RIRE2",
        superfamily="Ty3_Gypsy",
        strand="-",
        divergence=0.08,
        tsd=tsd_5,
        tsd_status="mismatch_1",
        tsd_3=tsd_3,
        tier=tier,
    )

    # Case 3: Degraded TSD (partial deletion)
    raw_tsd = make_tsd(5, seed=SEED_BASE + 403)
    tsd_5, tsd_3 = degrade_tsd(raw_tsd, mode="degraded", seed=SEED_BASE + 403)
    builder.insert_solo_ltr(
        position=5000,
        ltr_seq=ltr_rire2,
        family="RIRE2",
        superfamily="Ty3_Gypsy",
        strand="+",
        divergence=0.15,
        tsd=tsd_5,
        tsd_status="degraded",
        tsd_3=tsd_3,
        tier=tier,
    )

    # Case 4: Solo LTR with small internal remnant attached (50 bp fragment)
    #         Simulates incomplete deletion — common confounder
    remnant = generate_background(50, gc_content=0.42, seed=SEED_BASE + 404)
    ltr_with_remnant = ltr_rire2 + remnant   # partial internal region still attached
    tsd4 = make_tsd(5, seed=SEED_BASE + 404)
    builder.insert_solo_ltr(
        position=7500,
        ltr_seq=ltr_with_remnant,
        family="RIRE2",
        superfamily="Ty3_Gypsy",
        strand="+",
        divergence=0.12,
        tsd=tsd4,
        tsd_status="perfect",
        tsd_3=tsd4,
        tier=tier,
    )

    seq = builder.get_sequence()
    elements = builder.get_elements()

    write_fasta([("chr1 t0_edge_cases", seq)], outdir / "t0_edge_cases.fa")
    write_gff3(elements, outdir / "t0_edge_cases.gff3")
    write_bed(elements, outdir / "t0_edge_cases.bed")
    write_manifest(elements, outdir / "t0_edge_cases.json")

    print(f"  t0_edge_cases: {len(elements)} elements in {len(seq):,} bp genome")


def main():
    parser = argparse.ArgumentParser(description="Generate Tier 0 unit test sequences")
    parser.add_argument("--libdir", default="data/lib", help="Library directory")
    parser.add_argument("--outdir", default=OUTDIR_DEFAULT, help="Output directory")
    parser.add_argument("--tier", type=int, default=0, help="Tier label for GFF3")
    args = parser.parse_args()

    lib = Path(args.libdir)
    ltr_rire2 = _load_ltr(lib / "ltr_only.fa", "RIRE2_LTR")
    ltr_crm1 = _load_ltr(lib / "ltr_only.fa", "CRM1_LTR")
    internal_rire2 = _load_internal(lib / "internal_only.fa", "RIRE2")
    outdir = Path(args.outdir)

    print("Generating Tier 0 unit tests:")
    make_t0_clean(ltr_rire2, outdir / "clean", args.tier)
    make_t0_diverged(ltr_rire2, outdir / "diverged", args.tier)
    make_t0_intact_nearby(ltr_rire2, internal_rire2, outdir / "intact_nearby", args.tier)
    make_t0_edge_cases(ltr_rire2, ltr_crm1, outdir / "edge_cases", args.tier)
    print("Done.")


if __name__ == "__main__":
    main()
