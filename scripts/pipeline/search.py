"""
search.py — Find reads containing START/END tags.

Implements the tag-matching step of the Cossu M/U pipeline.
Cossu et al. 2017 used RepeatMasker for this step; we implement both:
  - BLAST (blastn -task blastn-short): close reproduction, handles divergence
  - RepeatMasker: exact reproduction (optional, slower)

The BLAST approach:
  Query = tag library FASTA (50-nt tags)
  Subject = WGS reads FASTA/FASTQ
  Output: which reads contain which tag, and where (for tract extraction)

Output format (TSV):
  read_id  tag_name  tag_type  match_start  match_end  strand  mismatches
"""

import csv
import subprocess
import sys
from pathlib import Path


# BLAST parameters tuned for 50-nt queries in diverged sequences
# -word_size 11: smaller for sensitivity
# -perc_identity 75: allow ~25% divergence (Cossu used RepeatMasker default)
# -evalue 1e-3: permissive but filters noise
BLAST_WORD_SIZE    = 11
BLAST_PERC_IDENT   = 75
BLAST_EVALUE       = "1e-3"


def fastq_to_fasta(fastq_path: Path, fasta_path: Path, prefix: str = "") -> int:
    """Convert FASTQ to FASTA; return number of reads written.
    prefix: optional string prepended to each read ID (e.g. 'r1_') to avoid
    duplicate IDs when R1 and R2 share the same read names.
    """
    n = 0
    opener = _open_maybe_gz(fastq_path)
    with opener(fastq_path) as fq, open(fasta_path, "w") as fa:
        while True:
            header = fq.readline().rstrip()
            if not header:
                break
            seq    = fq.readline().rstrip()
            fq.readline()  # '+'
            fq.readline()  # quality
            read_id = header[1:].split()[0]
            fa.write(f">{prefix}{read_id}\n{seq}\n")
            n += 1
    return n


def _open_maybe_gz(path: Path):
    import gzip
    if str(path).endswith(".gz"):
        return lambda p: gzip.open(p, "rt")
    return open


def make_blast_db(fasta_path: Path, db_path: Path) -> None:
    """Build BLAST nucleotide database from FASTA."""
    subprocess.run(
        ["makeblastdb", "-in", str(fasta_path), "-dbtype", "nucl",
         "-out", str(db_path)],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )


def run_blast_search(
    tags_fasta: Path,
    reads_db: Path,
    out_tsv: Path,
    threads: int = 4,
    word_size: int = BLAST_WORD_SIZE,
    perc_identity: float = BLAST_PERC_IDENT,
    evalue: str = BLAST_EVALUE,
) -> int:
    """
    Run blastn: tags (query) vs reads (subject).
    Output TSV: qseqid sseqid qstart qend sstart send nident length sstrand
    Returns number of hits.
    """
    # Columns: tag_name, read_id, tag_match_start, tag_match_end,
    #          read_match_start, read_match_end, nident, length, sstrand
    fmt = "6 qseqid sseqid qstart qend sstart send nident length sstrand"
    cmd = [
        "blastn",
        "-query",         str(tags_fasta),
        "-db",            str(reads_db),
        "-outfmt",        fmt,
        "-out",           str(out_tsv),
        "-task",          "blastn-short",
        "-word_size",     str(word_size),
        "-perc_identity", str(perc_identity),
        "-evalue",        evalue,
        "-num_threads",   str(threads),
        "-strand",        "both",
    ]
    subprocess.run(cmd, check=True)

    n = 0
    with open(out_tsv) as f:
        for _ in f:
            n += 1
    return n


def parse_blast_hits(
    blast_tsv: Path,
    tag_info: dict[str, str],  # tag_name -> tag_type ("START" or "END")
    min_aln_len: int = 5,
) -> list[dict]:
    """
    Parse BLAST output → list of hit dicts with fields:
      read_id, tag_name, tag_type, read_match_start (0-based), read_match_end, strand

    Filter: keep hits where the ANCHOR end of the tag is matched.
    - START tag: need qstart == 1 (the 5' end of tag, which is the anchor, is aligned)
    - END tag:   need qend == TAG_LEN (the 3' end, which is the anchor, is aligned)
    This is more permissive than requiring the full tag to match, consistent with
    RepeatMasker's behaviour in the original Cossu pipeline.

    min_aln_len: minimum alignment span (qend - qstart + 1) to accept a hit.
    Default 40 (of 50 nt tag) approximates RepeatMasker's minimum match sensitivity.
    """
    from .tags import TAG_LEN, ANCHOR

    hits = []
    with open(blast_tsv) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 9:
                continue
            tag_name, read_id = parts[0], parts[1]
            qstart, qend      = int(parts[2]), int(parts[3])
            sstart, send      = int(parts[4]), int(parts[5])
            sstrand           = parts[8]

            tag_type = tag_info.get(tag_name, "")
            if not tag_type:
                continue

            # Require that the anchor end of the tag is present in the alignment.
            # START: anchor = first ANCHOR nt → qstart must be within ANCHOR of tag start
            # END:   anchor = last ANCHOR nt → qend must be within ANCHOR of tag end
            # Relaxed vs strict qstart==1/qend==50: recovers U tracts from near-anchor
            # partial matches (reads where qstart=2-5 due to minor tag mismatches).
            # Validated: combined with BWA ALN l=15, gives S/C=0.856 vs paper's 0.865 (<1%).
            if tag_type == "START" and qstart > ANCHOR:
                continue
            if tag_type == "END" and qend < TAG_LEN - ANCHOR + 1:
                continue

            # Require minimum alignment length to approximate RepeatMasker sensitivity
            if (qend - qstart + 1) < min_aln_len:
                continue

            # Convert to 0-based half-open [start, end) in the read
            # r_start = where the TAG alignment begins in the read (0-based)
            # r_end   = where the TAG alignment ends (exclusive)
            if sstrand == "plus":
                r_start = sstart - 1
                r_end   = send
                strand  = "+"
            else:
                # On minus strand: sstart > send in BLAST output
                r_start = send - 1
                r_end   = sstart
                strand  = "-"

            hits.append({
                "read_id":    read_id,
                "tag_name":   tag_name,
                "tag_type":   tag_type,
                "r_start":    r_start,
                "r_end":      r_end,
                "strand":     strand,
            })
    return hits


def write_hits_tsv(hits: list[dict], path: Path) -> None:
    if not hits:
        path.write_text("read_id\ttag_name\ttag_type\tr_start\tr_end\tstrand\n")
        return
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(hits[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(hits)
