"""
Generate synthetic LTR-RT consensus sequences for simulation.

Real family sequences (RIRE2, CRM1, BARE-1, Athila6, ANGELA) are in Repbase
(licensed). We create synthetic representatives with correct structural properties:
  - TG...CA terminal dinucleotides (99% of plant LTR-RTs)
  - Realistic LTR lengths per family
  - PBS (primer binding site, 18 nt) immediately downstream of 5'LTR
  - PPT (polypurine tract, 15 nt) immediately upstream of 3'LTR
  - GC content 40-45% (plant LTR-RT typical)
  - Families are sufficiently diverged (≥40% pairwise) to allow family-specific detection
  - Each family has a fixed seed so sequences are reproducible

Outputs (data/lib/):
  ltr_families.fa    — full LTR-RT elements: 5'LTR + internal + 3'LTR
  ltr_only.fa        — just the LTR sequence (LTR pair identical at time 0)
  internal_only.fa   — just the internal region (gag-pol-int domain placeholder)
  background_tes.fa  — non-LTR TE fragments (LINE, SINE, DNA TE)
"""

import random
import argparse
from pathlib import Path

# ---------------------------------------------------------------------------
# Family definitions (name, superfamily, ltr_len, internal_len, gc, seed)
# LTR lengths from literature:
#   RIRE2 Ty3/Gypsy rice ~550 bp
#   CRM1  Ty3/Gypsy maize ~350 bp
#   BARE-1 Ty1/Copia barley ~1800 bp
#   Athila6 Ty3/Gypsy Arabidopsis ~530 bp
#   ANGELA Ty1/Copia wheat ~440 bp
# Internal region lengths (approximate, simplified):
#   Ty3/Gypsy: gag+pol+integrase ~5000-7000 bp
#   Ty1/Copia: gag+pol ~4000-6000 bp
# ---------------------------------------------------------------------------

FAMILIES = [
    {
        "name": "RIRE2",
        "superfamily": "Ty3_Gypsy",
        "classification": "LTR/Gypsy",
        "ltr_len": 550,
        "internal_len": 5800,
        "gc": 0.43,
        "seed": 1001,
    },
    {
        "name": "CRM1",
        "superfamily": "Ty3_Gypsy",
        "classification": "LTR/Gypsy",
        "ltr_len": 350,
        "internal_len": 5200,
        "gc": 0.41,
        "seed": 1002,
    },
    {
        "name": "BARE-1",
        "superfamily": "Ty1_Copia",
        "classification": "LTR/Copia",
        "ltr_len": 1800,
        "internal_len": 4500,
        "gc": 0.40,
        "seed": 1003,
    },
    {
        "name": "Athila6",
        "superfamily": "Ty3_Gypsy",
        "classification": "LTR/Gypsy",
        "ltr_len": 530,
        "internal_len": 6200,
        "gc": 0.44,
        "seed": 1004,
    },
    {
        "name": "ANGELA",
        "superfamily": "Ty1_Copia",
        "classification": "LTR/Copia",
        "ltr_len": 440,
        "internal_len": 4800,
        "gc": 0.42,
        "seed": 1005,
    },
]

# Non-LTR TE definitions for background noise (Tier 2+)
BACKGROUND_TE_FAMILIES = [
    {"name": "LINE1_synth", "classification": "LINE/L1", "length": 5000, "gc": 0.40, "seed": 2001},
    {"name": "LINE2_synth", "classification": "LINE/L1", "length": 3500, "gc": 0.38, "seed": 2002},
    {"name": "SINE1_synth", "classification": "SINE/tRNA", "length": 280, "gc": 0.42, "seed": 2003},
    {"name": "SINE2_synth", "classification": "SINE/tRNA", "length": 320, "gc": 0.44, "seed": 2004},
    {"name": "DNA_TE1_synth", "classification": "DNA/Tc1", "length": 800, "gc": 0.40, "seed": 2005},
    {"name": "DNA_TE2_synth", "classification": "DNA/hAT", "length": 1200, "gc": 0.43, "seed": 2006},
]


def gc_random_seq(length: int, gc: float, rng: random.Random) -> str:
    """Generate a random sequence with given GC content.

    The sequence is built from a pool of nucleotides with frequencies:
      G = gc/2, C = gc/2, A = (1-gc)/2, T = (1-gc)/2
    """
    at = 1.0 - gc
    pool = (
        ["G"] * int(gc / 2 * 1000)
        + ["C"] * int(gc / 2 * 1000)
        + ["A"] * int(at / 2 * 1000)
        + ["T"] * int(at / 2 * 1000)
    )
    return "".join(rng.choices(pool, k=length))


def make_ltr(ltr_len: int, gc: float, rng: random.Random) -> str:
    """Build a synthetic LTR with TG..CA termini and internal random sequence.

    Structure:
      TG + random_interior + CA
    where len(TG + interior + CA) = ltr_len.
    """
    interior_len = ltr_len - 4  # 2 nt at each end
    interior = gc_random_seq(interior_len, gc, rng)
    return "TG" + interior + "CA"


def make_pbs(rng: random.Random) -> str:
    """18-nt primer binding site (AT-rich, tRNA-complementary region placeholder)."""
    # PBS is AT-rich (~30% GC in real elements)
    return gc_random_seq(18, 0.30, rng)


def make_ppt(rng: random.Random) -> str:
    """15-nt polypurine tract (predominantly A/G)."""
    purines = ["A", "A", "A", "G", "G"]
    return "".join(rng.choices(purines, k=15))


def make_full_element(family: dict) -> tuple[str, str, str]:
    """Return (full_element, ltr_seq, internal_seq) for a family.

    Structure:
      5'LTR [TSD_placeholder] PBS internal PPT [TSD_placeholder] 3'LTR
    Returns:
      full_element  — complete LTR-RT (5'LTR + PBS + internal + PPT + 3'LTR)
      ltr_seq       — just the LTR sequence (same for 5' and 3' at time 0)
      internal_seq  — PBS + internal + PPT
    """
    rng = random.Random(family["seed"])
    ltr = make_ltr(family["ltr_len"], family["gc"], rng)
    pbs = make_pbs(rng)
    internal = gc_random_seq(family["internal_len"], family["gc"], rng)
    ppt = make_ppt(rng)
    full = ltr + pbs + internal + ppt + ltr
    internal_with_signals = pbs + internal + ppt
    return full, ltr, internal_with_signals


def write_fasta(records: list[tuple[str, str]], path: Path) -> None:
    """Write FASTA file. records = [(header, seq), ...]"""
    with open(path, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i : i + 60] + "\n")


def main():
    parser = argparse.ArgumentParser(description="Generate synthetic LTR-RT consensus library")
    parser.add_argument("--outdir", default="data/lib", help="Output directory")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    full_records = []
    ltr_records = []
    internal_records = []

    for fam in FAMILIES:
        full, ltr, internal = make_full_element(fam)
        name = fam["name"]
        cls = fam["classification"]

        # RepeatMasker-style header: >name#class/superfamily
        full_records.append(
            (
                f"{name}#LTR/{fam['superfamily']} [synthetic; LTR={fam['ltr_len']}bp internal={fam['internal_len']}bp]",
                full,
            )
        )
        ltr_records.append(
            (
                f"{name}_LTR#{cls} [synthetic; len={fam['ltr_len']}bp gc={fam['gc']}]",
                ltr,
            )
        )
        internal_records.append(
            (
                f"{name}_INT#{cls} [synthetic internal region; len={len(internal)}bp]",
                internal,
            )
        )

        print(
            f"  {name:10s}  LTR={len(ltr):4d} bp  internal={len(internal):5d} bp  "
            f"full={len(full):6d} bp  GC={fam['gc']:.0%}"
        )

    write_fasta(full_records, outdir / "ltr_families.fa")
    write_fasta(ltr_records, outdir / "ltr_only.fa")
    write_fasta(internal_records, outdir / "internal_only.fa")
    print(f"Wrote: {outdir}/ltr_families.fa, ltr_only.fa, internal_only.fa")

    # Background non-LTR TE fragments
    bg_records = []
    for bg in BACKGROUND_TE_FAMILIES:
        rng = random.Random(bg["seed"])
        seq = gc_random_seq(bg["length"], bg["gc"], rng)
        bg_records.append(
            (
                f"{bg['name']}#{bg['classification']} [synthetic background; len={bg['length']}bp]",
                seq,
            )
        )
        print(f"  {bg['name']:20s}  len={bg['length']:5d} bp  GC={bg['gc']:.0%}")

    write_fasta(bg_records, outdir / "background_tes.fa")
    print(f"Wrote: {outdir}/background_tes.fa")


if __name__ == "__main__":
    main()
