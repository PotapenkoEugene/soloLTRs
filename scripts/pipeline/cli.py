"""
soloLTRs pipeline CLI entry point.

Commands:
  sololtrs prepare  — build collapsed consensus + index
  sololtrs run      — full single-sample pipeline (align → depth → ratio)
  sololtrs batch    — run on many samples in parallel
  sololtrs pop      — aggregate per-sample results into population matrix
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from scripts.pipeline import annotate, align, depth, ratio, population
from scripts.pipeline.config import DEFAULT_THREADS


# ---------------------------------------------------------------------------
# prepare
# ---------------------------------------------------------------------------

def cmd_prepare(args):
    """Build collapsed consensus library and bwa-mem2 index."""
    out_fa = Path(args.out_dir) / "collapsed_consensus.fa"
    out_regions = Path(args.out_dir) / "regions.tsv"
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)

    if args.ltr and args.internal:
        regions = annotate.build_from_split(
            Path(args.ltr), Path(args.internal), out_fa, out_regions
        )
    else:
        regions = annotate.build_from_full(Path(args.lib), out_fa, out_regions)

    align.index_library(out_fa, force=args.reindex)
    print(f"\nLibrary prepared: {out_fa}")
    print(f"Regions map:      {out_regions}")


# ---------------------------------------------------------------------------
# run (single sample)
# ---------------------------------------------------------------------------

def run_sample(
    lib_dir: Path,
    r1: Path,
    r2: Path | None,
    out_dir: Path,
    threads: int,
    sample_id: str | None,
    bootstrap: bool,
    trim: int,
) -> Path:
    """Full pipeline for one sample. Returns path to ratios.tsv."""
    out_dir.mkdir(parents=True, exist_ok=True)
    consensus_fa = lib_dir / "collapsed_consensus.fa"
    regions_tsv = lib_dir / "regions.tsv"

    # Step 1: Align
    bam_path = out_dir / "aligned.bam"
    align.align_sample(
        consensus_fa=consensus_fa,
        r1=r1,
        r2=r2,
        out_bam=bam_path,
        threads=threads,
        sample_id=sample_id or out_dir.name,
    )

    # Step 2: Depth
    regions = annotate.load_regions(regions_tsv)
    depth_results = depth.compute_depth(bam_path, regions, trim=trim)
    depth_tsv = out_dir / "depth_table.tsv"
    depth.write_depth_table(depth_results, depth_tsv)

    # Step 3: Ratio
    ratio_results = []
    for d in depth_results:
        r = ratio.estimate_ratio(d)
        if bootstrap and d.d_ltr >= 10.0:
            r.ci_lo, r.ci_hi, r.sc_ci_lo, r.sc_ci_hi = ratio.bootstrap_ci(d)
        ratio_results.append(r)

    ratios_tsv = out_dir / "ratios.tsv"
    ratio.write_ratios(ratio_results, ratios_tsv)

    # Print summary
    n_detected = sum(1 for r in ratio_results if r.solo_detected)
    print(f"\n[{sample_id or out_dir.name}] Solo LTRs detected in {n_detected}/{len(ratio_results)} families")
    for r in ratio_results:
        if r.call not in ("absent", "low_coverage"):
            pav = "+" if r.solo_detected else "-"
            sc_str = f"{r.SC:.1f}" if r.SC == r.SC and r.SC != float("inf") else "inf"
            print(f"  [{pav}] {r.family:12s}  R={r.R:.3f}  S/C={sc_str:>6s}  ({r.confidence})")

    return ratios_tsv


def cmd_run(args):
    lib_dir = Path(args.lib_dir)
    r2 = Path(args.r2) if args.r2 else None
    run_sample(
        lib_dir=lib_dir,
        r1=Path(args.r1),
        r2=r2,
        out_dir=Path(args.out),
        threads=args.threads,
        sample_id=args.sample_id,
        bootstrap=args.bootstrap,
        trim=args.trim,
    )


# ---------------------------------------------------------------------------
# batch
# ---------------------------------------------------------------------------

def _run_one(kwargs):
    """Top-level function for multiprocessing (must be picklable)."""
    try:
        return run_sample(**kwargs), None
    except Exception as e:
        return None, str(e)


def cmd_batch(args):
    lib_dir = Path(args.lib_dir)
    out_root = Path(args.out)

    # Parse samples TSV
    samples = []
    with open(args.samples) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            samples.append(row)

    print(f"Batch: {len(samples)} samples, {args.parallel} parallel workers")

    jobs = []
    for row in samples:
        sid = row["sample_id"]
        r2 = Path(row["r2"]) if row.get("r2") else None
        jobs.append(dict(
            lib_dir=lib_dir,
            r1=Path(row["r1"]),
            r2=r2,
            out_dir=out_root / sid,
            threads=max(1, args.threads // args.parallel),
            sample_id=sid,
            bootstrap=args.bootstrap,
            trim=args.trim,
        ))

    failed = []
    with ProcessPoolExecutor(max_workers=args.parallel) as executor:
        futures = {executor.submit(_run_one, j): j["sample_id"] for j in jobs}
        for future in as_completed(futures):
            sid = futures[future]
            _, err = future.result()
            if err:
                print(f"  ERROR [{sid}]: {err}", file=sys.stderr)
                failed.append(sid)
            else:
                print(f"  DONE  [{sid}]")

    if failed:
        print(f"\nFailed samples: {', '.join(failed)}", file=sys.stderr)

    print(f"\nResults written to: {out_root}/")


# ---------------------------------------------------------------------------
# pop (population aggregation)
# ---------------------------------------------------------------------------

def cmd_pop(args):
    population.main.__module__  # trigger import check
    import sys as _sys
    _sys.argv = ["population"] + (
        ["--results-dir", args.results_dir] if args.results_dir else
        ["--samples", args.samples]
    ) + ["--out", args.out]
    population.main()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        prog="sololtrs",
        description="Reference-free solo LTR detection from short reads",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # --- prepare ---
    p_prep = sub.add_parser("prepare", help="Build collapsed consensus library")
    src = p_prep.add_mutually_exclusive_group(required=True)
    src.add_argument("--ltr", metavar="LTR_FA", help="ltr_only.fa")
    src.add_argument("--lib", metavar="FULL_FA", help="Full LTR-RT library FASTA")
    p_prep.add_argument("--internal", metavar="INT_FA", help="internal_only.fa (with --ltr)")
    p_prep.add_argument("--out-dir", required=True, help="Output directory for consensus + index")
    p_prep.add_argument("--reindex", action="store_true")
    p_prep.set_defaults(func=cmd_prepare)

    # --- run ---
    p_run = sub.add_parser("run", help="Run pipeline on a single sample")
    p_run.add_argument("--lib-dir", required=True, help="Library directory (from prepare)")
    p_run.add_argument("--r1", required=True, help="Forward reads FASTQ")
    p_run.add_argument("--r2", default=None, help="Reverse reads FASTQ")
    p_run.add_argument("--out", required=True, help="Output directory")
    p_run.add_argument("--threads", type=int, default=DEFAULT_THREADS)
    p_run.add_argument("--sample-id", default=None)
    p_run.add_argument("--bootstrap", action="store_true", help="Compute bootstrap CIs")
    p_run.add_argument("--trim", type=int, default=50, help="Edge trim bp (default: 50)")
    p_run.set_defaults(func=cmd_run)

    # --- batch ---
    p_bat = sub.add_parser("batch", help="Run pipeline on multiple samples")
    p_bat.add_argument("--lib-dir", required=True)
    p_bat.add_argument("--samples", required=True, help="TSV: sample_id, r1, r2")
    p_bat.add_argument("--out", required=True, help="Output root directory")
    p_bat.add_argument("--threads", type=int, default=DEFAULT_THREADS)
    p_bat.add_argument("--parallel", type=int, default=4, help="Samples in parallel")
    p_bat.add_argument("--bootstrap", action="store_true")
    p_bat.add_argument("--trim", type=int, default=50)
    p_bat.set_defaults(func=cmd_batch)

    # --- pop ---
    p_pop = sub.add_parser("pop", help="Aggregate samples into population matrix")
    src_pop = p_pop.add_mutually_exclusive_group(required=True)
    src_pop.add_argument("--results-dir", metavar="DIR")
    src_pop.add_argument("--samples", metavar="TSV")
    p_pop.add_argument("--out", required=True)
    p_pop.set_defaults(func=cmd_pop)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
