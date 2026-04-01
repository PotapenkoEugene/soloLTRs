"""
cli.py — Command-line interface for the Cossu M/U pipeline.

Usage:
  python -m scripts.pipeline.cli prepare --paralogs <FASTA> --out-dir <DIR>
  python -m scripts.pipeline.cli run     --reads <FASTQ[.gz]> --out-dir <DIR> [--threads N]
  python -m scripts.pipeline.cli all     --paralogs <FASTA> --reads <FASTQ[.gz]> --out-dir <DIR>
"""

import argparse
import csv
import sys
from pathlib import Path

from .tags import build_tag_fasta, load_fasta
from .search import (
    fastq_to_fasta, make_blast_db, run_blast_search, parse_blast_hits, write_hits_tsv
)
from .mu import (
    extract_tracts, build_bwa_index, map_tracts_bwa_aln,
    count_mu_from_bam, compute_sc_ratios, write_results_tsv,
)


def cmd_prepare(args):
    """
    Extract START/END tags from paralogs → write tags FASTA.
    Output: <out_dir>/tags.fa
    """
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paralogs_path = Path(args.paralogs)

    tags_fasta = out_dir / "tags.fa"
    tags = build_tag_fasta(paralogs_path, tags_fasta)
    print(f"Extracted {len(tags)} tags ({len(tags)//2} paralogs) → {tags_fasta}")

    # Save tag info TSV for downstream use
    tag_info_tsv = out_dir / "tag_info.tsv"
    with open(tag_info_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["tag_name", "tag_type"])
        for tag_name, tag_type, _ in tags:
            w.writerow([tag_name, tag_type])
    print(f"Tag info written → {tag_info_tsv}")


def cmd_run(args):
    """
    Run M/U pipeline:
      1. Convert reads FASTQ → FASTA (if needed)
      2. Build BLAST DB from reads FASTA
      3. BLAST: tags vs reads → hits TSV
      4. Parse hits → extract 20-nt tracts FASTA
      5. Build BWA index for paralogs
      6. BWA ALN: tracts → paralogs → BAM
      7. Count M/U from BAM
      8. Compute S/C ratios → results TSV
    """
    out_dir   = Path(args.out_dir)
    reads_inputs = [Path(r) for r in args.reads]
    threads   = args.threads

    # Locate prepare outputs
    tags_fasta   = out_dir / "tags.fa"
    tag_info_tsv = out_dir / "tag_info.tsv"
    paralogs_path = Path(args.paralogs)

    if not tags_fasta.exists() or not tag_info_tsv.exists():
        print(f"ERROR: Run 'prepare' first — {tags_fasta} or {tag_info_tsv} missing", file=sys.stderr)
        sys.exit(1)

    # Load tag info
    tag_info: dict[str, str] = {}
    with open(tag_info_tsv) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            tag_info[row["tag_name"]] = row["tag_type"]

    # Load paralogs (for BWA index + family listing)
    paralogs = load_fasta(paralogs_path)

    # ── Step 1: reads FASTQ → FASTA (all inputs concatenated) ────────────────
    reads_fasta = out_dir / "reads.fa"
    if not reads_fasta.exists():
        total_reads = 0
        first = True
        for i, reads_in in enumerate(reads_inputs):
            print(f"Converting {reads_in.name} ...", flush=True)
            tmp = out_dir / f"_tmp_{reads_in.stem}.fa"
            # Prefix read IDs with r1_, r2_, etc. to avoid duplicate IDs between
            # paired-end files (both R1 and R2 share the same read names).
            prefix = f"r{i+1}_" if len(reads_inputs) > 1 else ""
            n = fastq_to_fasta(reads_in, tmp, prefix=prefix)
            total_reads += n
            print(f"  {n:,} reads")
            mode = "w" if first else "a"
            with open(tmp) as src, open(reads_fasta, mode) as dst:
                for line in src:
                    dst.write(line)
            tmp.unlink()
            first = False
        print(f"  {total_reads:,} reads total → {reads_fasta}")
    else:
        print(f"Reads FASTA already exists: {reads_fasta}")

    # ── Step 2: BLAST DB from reads ──────────────────────────────────────────
    blast_db = out_dir / "reads_db" / "reads"
    blast_db.parent.mkdir(parents=True, exist_ok=True)
    db_marker = out_dir / "reads_db" / "reads.nin"
    if not db_marker.exists():
        print(f"Building BLAST DB from reads ...", flush=True)
        make_blast_db(reads_fasta, blast_db)
        print("  Done")
    else:
        print(f"BLAST DB already exists: {blast_db}")

    # ── Step 3: BLAST search ──────────────────────────────────────────────────
    blast_tsv = out_dir / "blast_hits.tsv"
    if not blast_tsv.exists():
        print(f"Running BLAST: tags vs reads ...", flush=True)
        n_hits = run_blast_search(tags_fasta, blast_db, blast_tsv, threads=threads)
        print(f"  {n_hits:,} raw BLAST hits")
    else:
        print(f"BLAST output already exists: {blast_tsv}")

    # ── Step 4: Parse hits ────────────────────────────────────────────────────
    hits_tsv = out_dir / "hits.tsv"
    if not hits_tsv.exists():
        print("Parsing BLAST hits ...", flush=True)
        hits = parse_blast_hits(blast_tsv, tag_info)
        write_hits_tsv(hits, hits_tsv)
        print(f"  {len(hits):,} valid hits → {hits_tsv}")
    else:
        print(f"Hits TSV already exists: {hits_tsv}")

    # ── Step 5: Extract 20-nt tracts ─────────────────────────────────────────
    tracts_fasta = out_dir / "tracts.fa"
    if not tracts_fasta.exists():
        print("Extracting 20-nt tracts ...", flush=True)
        n_tracts = extract_tracts(hits_tsv, reads_fasta, tracts_fasta)
        print(f"  {n_tracts:,} tracts")
    else:
        print(f"Tracts FASTA already exists: {tracts_fasta}")

    # ── Step 6: BWA index for paralogs ───────────────────────────────────────
    bwa_index_marker = paralogs_path.with_suffix(".fa.bwt")
    if not bwa_index_marker.exists():
        # BWA writes index files alongside the FASTA
        alt_marker = Path(str(paralogs_path) + ".bwt")
        if not alt_marker.exists():
            print(f"Building BWA index for paralogs ...", flush=True)
            build_bwa_index(paralogs_path)
            print("  Done")
    else:
        print(f"BWA index already exists for {paralogs_path}")

    # ── Step 7: BWA ALN mapping ───────────────────────────────────────────────
    bam_path = out_dir / "tracts.bam"
    if not bam_path.exists():
        print("Mapping tracts with BWA ALN ...", flush=True)
        map_tracts_bwa_aln(tracts_fasta, paralogs_path, bam_path, threads=threads)
        print(f"  BAM written → {bam_path}")
    else:
        print(f"BAM already exists: {bam_path}")

    # ── Step 8: Count M/U + compute S/C ──────────────────────────────────────
    print("Counting M/U ...", flush=True)
    counts = count_mu_from_bam(bam_path, hits_tsv, paralogs)
    rows = compute_sc_ratios(counts)

    results_tsv = out_dir / "results.tsv"
    write_results_tsv(rows, results_tsv)
    print(f"\nResults → {results_tsv}\n")

    # Print summary table
    header = f"{'Family':<30} {'M':>8} {'U':>8} {'M+U':>8} {'U/M':>8} {'S/C':>8}"
    print(header)
    print("-" * len(header))
    for row in rows:
        print(
            f"{row['family']:<30} {row['M']:>8} {row['U']:>8} "
            f"{row['M_plus_U']:>8} {row['U_over_M']:>8} {row['S_to_C']:>8}"
        )


def cmd_all(args):
    cmd_prepare(args)
    cmd_run(args)


def main():
    parser = argparse.ArgumentParser(
        description="Cossu M/U solo:intact ratio pipeline",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # ── prepare ───────────────────────────────────────────────────────────────
    p_prep = sub.add_parser("prepare", help="Extract START/END tags from paralogs")
    p_prep.add_argument("--paralogs", required=True, help="Paralog library FASTA")
    p_prep.add_argument("--out-dir",  required=True, help="Output directory")

    # ── run ───────────────────────────────────────────────────────────────────
    p_run = sub.add_parser("run", help="Run M/U pipeline (needs prepare output)")
    p_run.add_argument("--paralogs", required=True, help="Paralog library FASTA")
    p_run.add_argument("--reads",    required=True, nargs="+",
                       help="WGS reads FASTQ or FASTQ.gz (R1 [R2 ...])")
    p_run.add_argument("--out-dir",  required=True, help="Output directory")
    p_run.add_argument("--threads",  type=int, default=4, help="CPU threads (default: 4)")

    # ── all ───────────────────────────────────────────────────────────────────
    p_all = sub.add_parser("all", help="Run prepare + run in one step")
    p_all.add_argument("--paralogs", required=True, help="Paralog library FASTA")
    p_all.add_argument("--reads",    required=True, nargs="+",
                       help="WGS reads FASTQ or FASTQ.gz (R1 [R2 ...])")
    p_all.add_argument("--out-dir",  required=True, help="Output directory")
    p_all.add_argument("--threads",  type=int, default=4, help="CPU threads (default: 4)")

    args = parser.parse_args()
    if args.command == "prepare":
        cmd_prepare(args)
    elif args.command == "run":
        cmd_run(args)
    elif args.command == "all":
        cmd_all(args)


if __name__ == "__main__":
    main()
