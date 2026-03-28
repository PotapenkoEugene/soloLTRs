"""
Tier 2: Divergence × background matrix (~50 Mb, 3 chromosomes).

Purpose: Sensitivity-vs-divergence curve; specificity against non-LTR TEs.

Genome specs:
  - 50 Mb, 3 chromosomes (10 + 20 + 20 Mb)
  - GC composition:
      chr1: 42% uniform
      chr2: 42% base, 3 AT-rich islands (30% GC, 100 kb each)
      chr3: 42% base, 2 GC-rich islands (55% GC, 100 kb each)
  - 200 solo LTRs from 5 families (RIRE2, CRM1, BARE-1, Athila6, ANGELA):
      50 at 0-5% divergence
      50 at 5-10% divergence
      50 at 10-20% divergence
      50 at 20-30% divergence
      (10 per family per divergence bin)
  - 30 intact LTR-RTs (6 per family)
  - 100 non-LTR background fragments (from background_tes.fa)
  - TSDs: 70% perfect, 20% mismatch_1, 10% degraded

Outputs:
  genome/t2_genome.fa              — gitignored
  ground_truth/t2_ground_truth.gff3
  ground_truth/t2_ground_truth.bed
  ground_truth/t2_ground_truth.json
"""

import argparse
import random
import sys
from pathlib import Path
from collections import Counter

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


OUTDIR_DEFAULT = "data/sim/t2"
SEED = 200

# Chromosome sizes
CHR_SIZES = {"chr1": 10_000_000, "chr2": 20_000_000, "chr3": 20_000_000}

FAMILIES = ["RIRE2", "CRM1", "BARE-1", "Athila6", "ANGELA"]

# Divergence bins: (label, min, max)
DIV_BINS = [
    ("0-5", 0.00, 0.05),
    ("5-10", 0.05, 0.10),
    ("10-20", 0.10, 0.20),
    ("20-30", 0.20, 0.30),
]
SOLOS_PER_FAMILY_PER_BIN = 10   # 5 families × 4 bins × 10 = 200 total
INTACT_PER_FAMILY = 6           # 5 families × 6 = 30 total
N_BACKGROUND_TES = 100

# TSD mode probabilities: 70% perfect, 20% mismatch_1, 10% degraded
TSD_MODES = ["perfect"] * 70 + ["mismatch_1"] * 20 + ["degraded"] * 10


def _load_seqs(path: Path, tag: str = "") -> dict[str, str]:
    """Load all sequences from a FASTA. Key = first word after '>' stripped of tag."""
    seqs: dict[str, list[str]] = {}
    current = None
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                raw = line[1:].split("#")[0].split()[0]
                current = raw.replace(tag, "").strip("_")
                seqs[current] = []
            elif current:
                seqs[current].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def make_gc_patched_background(
    total_len: int,
    base_gc: float,
    patches: list[tuple[int, int, float]],  # (start, end, gc)
    seed: int,
) -> str:
    """Generate background sequence with GC-composition patches.

    patches: list of (start_pos, end_pos, gc_content) — positions are
    within the final sequence (0-based, non-overlapping, sorted).
    Regions outside patches get base_gc.
    """
    seq_list = list(generate_background(total_len, base_gc, seed))
    for start, end, gc in patches:
        patch_len = end - start
        patch = generate_background(patch_len, gc, seed=seed + start)
        seq_list[start:end] = list(patch)
    return "".join(seq_list)


def place_elements_on_chromosome(
    genome_len: int,
    n_elements: int,
    min_spacing: int,
    rng: random.Random,
) -> list[int]:
    """Return n_elements positions spaced at least min_spacing apart."""
    slots = list(range(min_spacing, genome_len - min_spacing, min_spacing))
    if len(slots) < n_elements:
        raise ValueError(
            f"Cannot fit {n_elements} elements with min_spacing={min_spacing} in {genome_len} bp"
        )
    chosen = sorted(rng.sample(slots, n_elements))
    return chosen


def make_t2(lib: Path, outdir: Path) -> None:
    rng = random.Random(SEED)

    # Load LTR and internal sequences
    ltr_seqs = _load_seqs(lib / "ltr_only.fa", "_LTR")
    int_seqs = _load_seqs(lib / "internal_only.fa", "_INT")
    bg_seqs = _load_seqs(lib / "background_tes.fa")

    # --- Build each chromosome ---
    chr_seqs: dict[str, str] = {}
    chr_builders: dict[str, GenomeBuilder] = {}

    # chr1: uniform 42% GC
    chr_seqs["chr1"] = generate_background(CHR_SIZES["chr1"], 0.42, seed=SEED + 1)
    chr_builders["chr1"] = GenomeBuilder(chr_seqs["chr1"], chrom="chr1", tier=2)

    # chr2: 3 AT-rich islands at ~3Mb, ~10Mb, ~16Mb
    chr_seqs["chr2"] = make_gc_patched_background(
        CHR_SIZES["chr2"], 0.42,
        patches=[(3_000_000, 3_100_000, 0.30),
                 (10_000_000, 10_100_000, 0.30),
                 (16_000_000, 16_100_000, 0.30)],
        seed=SEED + 2,
    )
    chr_builders["chr2"] = GenomeBuilder(chr_seqs["chr2"], chrom="chr2", tier=2)

    # chr3: 2 GC-rich islands at ~5Mb, ~14Mb
    chr_seqs["chr3"] = make_gc_patched_background(
        CHR_SIZES["chr3"], 0.42,
        patches=[(5_000_000, 5_100_000, 0.55),
                 (14_000_000, 14_100_000, 0.55)],
        seed=SEED + 3,
    )
    chr_builders["chr3"] = GenomeBuilder(chr_seqs["chr3"], chrom="chr3", tier=2)

    chroms = list(chr_builders.keys())

    # --- Distribute elements across chromosomes ---
    # Assign each element a chromosome randomly (weighted by size)
    chr_weights = [CHR_SIZES[c] for c in chroms]
    total_size = sum(chr_weights)
    chr_probs = [w / total_size for w in chr_weights]

    def pick_chrom() -> str:
        return rng.choices(chroms, weights=chr_probs, k=1)[0]

    # Track positions already used per chromosome (to avoid overlap)
    used_positions: dict[str, list[tuple[int, int]]] = {c: [] for c in chroms}
    MAX_LTR_LEN = 1800   # BARE-1 LTR
    INTACT_MAX = 8200    # longest full element (BARE-1)
    MIN_SPACING = max(MAX_LTR_LEN, INTACT_MAX) + 500

    def find_position(chrom: str, element_len: int, attempts: int = 1000) -> int | None:
        chr_len = CHR_SIZES[chrom]
        margin = MIN_SPACING
        for _ in range(attempts):
            pos = rng.randint(margin, chr_len - element_len - margin)
            if not any(abs(pos - u[0]) < MIN_SPACING for u in used_positions[chrom]):
                used_positions[chrom].append((pos, pos + element_len))
                return pos
        return None

    def pick_tsd_mode(seed_val: int) -> tuple[str, str, str, str]:
        """Return (tsd_5, tsd_3, tsd_status, mode_used)."""
        mode = rng.choice(TSD_MODES)
        tsd_raw = make_tsd(5, seed=seed_val)
        tsd_5, tsd_3 = degrade_tsd(tsd_raw, mode=mode, seed=seed_val)
        return tsd_5, tsd_3, mode, mode

    # --- Insert 200 solo LTRs (5 families × 4 bins × 10) ---
    elem_seed = SEED + 10000
    solo_count = 0
    failed = 0
    for bin_label, div_min, div_max in DIV_BINS:
        for fam in FAMILIES:
            ltr_key = f"{fam}_LTR" if f"{fam}_LTR" in ltr_seqs else fam
            ltr_consensus = ltr_seqs.get(fam, ltr_seqs.get(f"{fam}_LTR", ""))
            if not ltr_consensus:
                print(f"WARNING: LTR not found for family {fam}")
                continue

            for _ in range(SOLOS_PER_FAMILY_PER_BIN):
                div = rng.uniform(div_min, div_max)
                ltr_mut = mutate_sequence(ltr_consensus, div, seed=elem_seed)
                chrom = pick_chrom()
                pos = find_position(chrom, len(ltr_mut) + 10)
                if pos is None:
                    failed += 1
                    elem_seed += 1
                    continue

                tsd_5, tsd_3, tsd_status, _ = pick_tsd_mode(elem_seed)
                strand = rng.choice(["+", "-"])
                chr_builders[chrom].insert_solo_ltr(
                    position=pos,
                    ltr_seq=ltr_mut,
                    family=fam,
                    superfamily="Ty3_Gypsy" if fam in ("RIRE2", "CRM1", "Athila6") else "Ty1_Copia",
                    strand=strand,
                    divergence=div,
                    divergence_bin=bin_label,
                    tsd=tsd_5,
                    tsd_status=tsd_status,
                    tsd_3=tsd_3,
                    tier=2,
                )
                solo_count += 1
                elem_seed += 1

    # --- Insert 30 intact LTR-RTs (5 families × 6) ---
    intact_count = 0
    for fam in FAMILIES:
        ltr_consensus = ltr_seqs.get(fam, "")
        internal_consensus = int_seqs.get(fam, "")
        if not ltr_consensus or not internal_consensus:
            continue
        full_len = 2 * len(ltr_consensus) + len(internal_consensus) + 10  # +TSD
        for _ in range(INTACT_PER_FAMILY):
            chrom = pick_chrom()
            pos = find_position(chrom, full_len)
            if pos is None:
                continue
            tsd_5, tsd_3, tsd_status, _ = pick_tsd_mode(elem_seed)
            strand = rng.choice(["+", "-"])
            div = rng.uniform(0.01, 0.05)
            chr_builders[chrom].insert_intact_ltr_rt(
                position=pos,
                ltr_seq=ltr_consensus,
                internal_seq=internal_consensus,
                family=fam,
                superfamily="Ty3_Gypsy" if fam in ("RIRE2", "CRM1", "Athila6") else "Ty1_Copia",
                strand=strand,
                divergence=div,
                tsd=tsd_5,
                tsd_status=tsd_status,
                tsd_3=tsd_3,
            )
            intact_count += 1
            elem_seed += 1

    # --- Insert 100 background non-LTR TE fragments ---
    bg_keys = list(bg_seqs.keys())
    bg_count = 0
    for _ in range(N_BACKGROUND_TES):
        bg_key = rng.choice(bg_keys)
        bg_seq_full = bg_seqs[bg_key]
        # Take a random fragment (200–2000 bp)
        frag_len = rng.randint(200, min(2000, len(bg_seq_full)))
        start_in_bg = rng.randint(0, len(bg_seq_full) - frag_len)
        frag = bg_seq_full[start_in_bg : start_in_bg + frag_len]

        chrom = pick_chrom()
        pos = find_position(chrom, frag_len)
        if pos is None:
            continue
        strand = rng.choice(["+", "-"])
        chr_builders[chrom].insert_background_te(
            position=pos,
            te_seq=frag,
            family=bg_key,
            strand=strand,
        )
        bg_count += 1
        elem_seed += 1

    # --- Collect all sequences and elements ---
    all_seqs: list[tuple[str, str]] = []
    all_elements = []
    for chrom in chroms:
        seq = chr_builders[chrom].get_sequence()
        all_seqs.append((f"{chrom} t2_genome", seq))
        all_elements.extend(chr_builders[chrom].get_elements())

    # --- Write outputs ---
    genome_dir = outdir / "genome"
    genome_dir.mkdir(parents=True, exist_ok=True)
    gt_dir = outdir / "ground_truth"
    gt_dir.mkdir(parents=True, exist_ok=True)

    write_fasta(all_seqs, genome_dir / "t2_genome.fa")
    write_gff3(all_elements, gt_dir / "t2_ground_truth.gff3")
    write_bed(all_elements, gt_dir / "t2_ground_truth.bed")
    write_manifest(all_elements, gt_dir / "t2_ground_truth.json")

    # --- Report ---
    solos = [e for e in all_elements if e.elem_type == "solo_LTR"]
    intact = [e for e in all_elements if e.elem_type == "intact_LTR_RT"]
    bg = [e for e in all_elements if e.elem_type == "background_TE"]

    total_len = sum(CHR_SIZES.values())
    print(f"Tier 2 genome: {sum(len(s) for _, s in all_seqs):,} bp ({len(chroms)} chromosomes)")
    print(f"  Solo LTRs:      {len(solos)} (target 200, failed placements: {failed})")
    print(f"  Intact LTR-RTs: {len(intact)} (target 30)")
    print(f"  Background TEs: {len(bg)} (target 100)")
    print(f"  Genome FASTA:   {genome_dir}/t2_genome.fa")
    print(f"  Ground truth:   {gt_dir}/")

    # Divergence bin distribution
    bin_counts = Counter(e.divergence_bin for e in solos)
    print("\n  Solo LTR counts by divergence bin:")
    for b in ["0-5", "5-10", "10-20", "20-30"]:
        print(f"    {b}: {bin_counts.get(b, 0)}")

    # Family distribution
    fam_counts = Counter(e.family for e in solos)
    print("\n  Solo LTR counts by family:")
    for fam in FAMILIES:
        print(f"    {fam}: {fam_counts.get(fam, 0)}")

    # TSD status distribution
    tsd_counts = Counter(e.tsd_status for e in solos)
    print("\n  TSD status distribution:")
    for status in ("perfect", "mismatch_1", "degraded"):
        print(f"    {status}: {tsd_counts.get(status, 0)}")


def main():
    parser = argparse.ArgumentParser(description="Generate Tier 2 divergence × background matrix")
    parser.add_argument("--libdir", default="data/lib")
    parser.add_argument("--outdir", default=OUTDIR_DEFAULT)
    args = parser.parse_args()
    make_t2(Path(args.libdir), Path(args.outdir))


if __name__ == "__main__":
    main()
