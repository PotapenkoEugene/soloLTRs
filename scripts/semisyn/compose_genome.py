"""
compose_genome.py — Build a virtual genome FASTA per sample.

For each sample in the manifest, constructs a small virtual genome by:
  1. Collecting the assigned loci sequences (with flanking) from loci_with_flanking.fa
  2. Concatenating them on a single virtual chromosome with 2 kb spacers
  3. Adding 50 kb of random background padding at start and end (for non-mapping reads)

Virtual genomes are intentionally small (typically 2-10 Mb per sample) because the
pipeline maps reads to the TE consensus library, not to the genome. The genome only
needs to contain the TE copies + flanking to generate realistic reads.

Also writes per-sample ground truth JSON (n_solo, n_intact, true S/C per family)
for use by evaluate_semisyn.py.

Usage:
    python scripts/semisyn/compose_genome.py \\
        --manifest data/semisyn/samples/t2_mock/manifest.json \\
        --loci-fa data/semisyn/extracted/t2_mock/loci_with_flanking.fa \\
        --loci-meta data/semisyn/extracted/t2_mock/loci_metadata.json \\
        --out-dir data/semisyn/samples/t2_mock
"""

import argparse
import json
import random
import sys
from pathlib import Path


BG_SPACER_LEN = 2000   # bp between loci
BG_PAD_LEN = 50_000    # bp at start and end of virtual chromosome
BG_SEED_BASE = 9999    # per-sample seed offset for background


def load_fasta(path: Path) -> dict[str, str]:
    """Load FASTA into {header_prefix -> sequence} using the first token of header."""
    seqs = {}
    current_id = None
    current_parts = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    seqs[current_id] = "".join(current_parts)
                current_id = line[1:].split()[0]  # first token after >
                current_parts = []
            else:
                current_parts.append(line.upper())
    if current_id is not None:
        seqs[current_id] = "".join(current_parts)
    return seqs


def generate_background(length: int, gc_content: float, seed: int) -> str:
    """Random DNA with controlled GC content."""
    rng = random.Random(seed)
    gc = gc_content / 2
    at = (1 - gc_content) / 2
    weights = [at, gc, gc, at]  # A C G T
    bases = "ACGT"
    return "".join(rng.choices(bases, weights=weights, k=length))


def write_fasta_record(f, header: str, seq: str, wrap: int = 60) -> None:
    f.write(f">{header}\n")
    for i in range(0, len(seq), wrap):
        f.write(seq[i:i+wrap] + "\n")


def compose_sample(
    sample: dict,
    loci_seqs: dict[str, str],
    bg_gc: float,
    sample_dir: Path,
    sample_seed: int,
) -> None:
    sample_id = sample["sample_id"]
    loci_ids = sample["loci"]

    # Build the virtual chromosome sequence
    chrom_parts = []

    # 5' background padding
    chrom_parts.append(generate_background(BG_PAD_LEN, bg_gc, sample_seed))

    # Loci with spacers
    missing = []
    for locus_id in loci_ids:
        seq = loci_seqs.get(locus_id)
        if seq is None:
            missing.append(locus_id)
            continue
        chrom_parts.append(seq)
        # Add spacer between loci (not after the last one before padding)
        chrom_parts.append(generate_background(BG_SPACER_LEN, bg_gc, sample_seed + hash(locus_id) % 10_000))

    if missing:
        print(f"  Warning: {len(missing)} loci missing from loci_with_flanking.fa: {missing[:5]}",
              file=sys.stderr)

    # 3' background padding
    chrom_parts.append(generate_background(BG_PAD_LEN, bg_gc, sample_seed + 1))

    chrom_seq = "".join(chrom_parts)

    # Write genome FASTA
    genome_path = sample_dir / "genome.fa"
    with open(genome_path, "w") as f:
        write_fasta_record(f, f"virtual_chr1 sample={sample_id}", chrom_seq)

    # Write per-sample ground truth JSON
    truth = {
        "sample_id": sample_id,
        "species": None,  # filled in by caller
        "truth_SC": sample["truth_SC"],
        "n_solo": sample["n_solo"],
        "n_intact": sample["n_intact"],
        "target_SC": sample["target_SC"],
        "n_loci": len(loci_ids),
        "genome_len": len(chrom_seq),
    }
    (sample_dir / "truth.json").write_text(json.dumps(truth, indent=2))

    return len(chrom_seq)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--loci-fa", required=True, type=Path,
                        help="loci_with_flanking.fa from extract_sequences.py")
    parser.add_argument("--loci-meta", required=True, type=Path,
                        help="loci_metadata.json from extract_sequences.py")
    parser.add_argument("--bg-gc", type=float, default=0.43,
                        help="Background GC content (default: 0.43 = rice-like)")
    parser.add_argument("--out-dir", required=True, type=Path,
                        help="Output directory (sample subdirs created here)")
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text())
    species = manifest.get("species", "unknown")

    print(f"Loading {args.loci_fa} ...", flush=True)
    loci_seqs = load_fasta(args.loci_fa)
    print(f"  Loaded {len(loci_seqs)} locus sequences")

    total_loci_used = 0
    for s_idx, sample in enumerate(manifest["samples"]):
        sample_id = sample["sample_id"]
        sample_dir = args.out_dir / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)

        n_loci = len(sample["loci"])
        genome_len = compose_sample(
            sample=sample,
            loci_seqs=loci_seqs,
            bg_gc=args.bg_gc,
            sample_dir=sample_dir,
            sample_seed=BG_SEED_BASE + s_idx * 100,
        )
        total_loci_used += n_loci

        # Fill in species in truth.json
        truth = json.loads((sample_dir / "truth.json").read_text())
        truth["species"] = species
        (sample_dir / "truth.json").write_text(json.dumps(truth, indent=2))

        families_summary = "  ".join(
            f"{f}:{sample['n_solo'].get(f,0)}/{sample['n_intact'].get(f,0)}"
            for f in manifest["families"] if f in sample["truth_SC"]
        )
        print(f"  {sample_id}: {n_loci} loci, genome {genome_len/1e6:.2f} Mb  [{families_summary}]")

    print(f"\nComposed {len(manifest['samples'])} virtual genomes in {args.out_dir}/")
    print(f"Total loci embedded: {total_loci_used}")


if __name__ == "__main__":
    main()
