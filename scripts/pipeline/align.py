"""
Align short reads to the collapsed consensus library using bwa-mem2.

Only TE-derived reads (a tiny fraction of total reads from a resequenced
genome) will map. The output BAM is therefore very small, making this step
fast even for 5 Gb genomes at 10x coverage.

Commands run:
  bwa-mem2 index collapsed_consensus.fa   (once per library)
  bwa-mem2 mem -t N collapsed_consensus.fa R1.fq.gz R2.fq.gz
    | samtools view -b -F 4 -q MIN_MAPQ
    | samtools sort -o sample.bam
  samtools index sample.bam
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

from scripts.pipeline.config import (
    DEFAULT_THREADS,
    MIN_MAPQ_ALIGN,
    BWA_MISMATCH_PENALTY,
    BWA_MIN_SCORE,
    BWA_MIN_SEED_LEN,
)


def check_tool(name: str) -> None:
    if shutil.which(name) is None:
        sys.exit(
            f"ERROR: '{name}' not found in PATH.\n"
            f"  Install via: nix shell nixpkgs#{name.replace('-', '')}"
        )


def index_library(consensus_fa: Path, force: bool = False) -> None:
    """Build bwa-mem2 index if not already present."""
    index_marker = consensus_fa.with_suffix(consensus_fa.suffix + ".bwt.2bit.64")
    if index_marker.exists() and not force:
        return
    check_tool("bwa-mem2")
    cmd = ["bwa-mem2", "index", str(consensus_fa)]
    print(f"[align] Indexing: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def align_sample(
    consensus_fa: Path,
    r1: Path,
    r2: Path | None,
    out_bam: Path,
    threads: int = DEFAULT_THREADS,
    min_mapq: int = MIN_MAPQ_ALIGN,
    sample_id: str | None = None,
) -> None:
    """Align reads to collapsed consensus, output sorted BAM (mapped only)."""
    check_tool("bwa-mem2")
    check_tool("samtools")

    out_bam.parent.mkdir(parents=True, exist_ok=True)

    # Build read-group tag if sample_id provided
    rg = ""
    if sample_id:
        rg = f"-R @RG\\tID:{sample_id}\\tSM:{sample_id}\\tPL:ILLUMINA"

    bwa_cmd = [
        "bwa-mem2", "mem",
        "-t", str(threads),
        "-k", str(BWA_MIN_SEED_LEN),
        "-B", str(BWA_MISMATCH_PENALTY),
        "-T", str(BWA_MIN_SCORE),
    ]
    if rg:
        bwa_cmd += ["-R", f"@RG\tID:{sample_id}\tSM:{sample_id}\tPL:ILLUMINA"]
    bwa_cmd += [str(consensus_fa), str(r1)]
    if r2 is not None:
        bwa_cmd.append(str(r2))

    view_cmd = [
        "samtools", "view",
        "-b",
        "-F", "4",           # discard unmapped reads
        "-q", str(min_mapq), # minimum MAPQ
        "-@", str(max(1, threads // 2)),
    ]

    sort_cmd = [
        "samtools", "sort",
        "-@", str(max(1, threads // 2)),
        "-o", str(out_bam),
    ]

    print(f"[align] Aligning {r1.name} → {out_bam.name}")
    print(f"[align] threads={threads}  min_mapq={min_mapq}")

    # Chain: bwa-mem2 | samtools view | samtools sort
    p1 = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE, stderr=sys.stderr)
    p2 = subprocess.Popen(view_cmd, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=sys.stderr)
    p1.stdout.close()
    p3 = subprocess.Popen(sort_cmd, stdin=p2.stdout, stderr=sys.stderr)
    p2.stdout.close()

    p1.wait()
    p2.wait()
    p3.wait()

    if p1.returncode != 0:
        sys.exit(f"ERROR: bwa-mem2 failed (exit {p1.returncode})")
    if p2.returncode != 0:
        sys.exit(f"ERROR: samtools view failed (exit {p2.returncode})")
    if p3.returncode != 0:
        sys.exit(f"ERROR: samtools sort failed (exit {p3.returncode})")

    # Index the BAM
    idx_cmd = ["samtools", "index", str(out_bam)]
    subprocess.run(idx_cmd, check=True)

    # Report mapped read count
    count_cmd = ["samtools", "view", "-c", str(out_bam)]
    result = subprocess.run(count_cmd, capture_output=True, text=True)
    n_mapped = result.stdout.strip()
    print(f"[align] Mapped reads: {n_mapped}")


def main():
    parser = argparse.ArgumentParser(description="Align reads to collapsed consensus library")
    parser.add_argument("--lib", required=True, help="Collapsed consensus FASTA")
    parser.add_argument("--r1", required=True, help="Forward reads FASTQ (can be .gz)")
    parser.add_argument("--r2", default=None, help="Reverse reads FASTQ (optional for PE)")
    parser.add_argument("--out", required=True, help="Output BAM path")
    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS)
    parser.add_argument("--sample-id", default=None, help="Sample ID for read group")
    parser.add_argument("--reindex", action="store_true", help="Force re-index of library")
    args = parser.parse_args()

    consensus_fa = Path(args.lib)
    r2 = Path(args.r2) if args.r2 else None

    index_library(consensus_fa, force=args.reindex)
    align_sample(
        consensus_fa=consensus_fa,
        r1=Path(args.r1),
        r2=r2,
        out_bam=Path(args.out),
        threads=args.threads,
        sample_id=args.sample_id,
    )


if __name__ == "__main__":
    main()
