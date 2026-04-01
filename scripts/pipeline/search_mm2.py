"""
search_mm2.py — Find reads containing START/END tags using minimap2.

Replaces the BLAST-based tag search in search.py.
Strategy: index the ~1KB tag library as a tiny reference and map WGS reads to it.
This is the reverse of BLAST (tags-as-query / reads-as-db), but produces identical
hits.tsv + tracts.fa so downstream BWA ALN + M/U counting is unchanged.

Key advantage: minimap2 reads FASTQ.gz directly, no FASTQ→FASTA conversion,
no makeblastdb, and tag-indexed search is far faster than BLAST on a large reads DB.
"""

import re
import subprocess
from pathlib import Path

from .search import write_hits_tsv
from .tags import ANCHOR, FLANK, TAG_LEN, extract_tract_from_read

# minimap2 parameters tuned for ~75% identity matches on 50-nt tags.
# A 50-nt alignment at 75% identity has 37 matches and 13 mismatches.
# With A=2, B=3: score = 37*2 - 13*3 = 35, well above MIN_SCORE=15.
# At 70% identity: score = 35*2 - 15*3 = 25, still passes → sensitive.
# Note: B=1 triggers minimap2's internal constraint (score-N must be < B),
# so B>=2 is required. B=3 (vs default B=4) increases sensitivity for diverged tags.
MM2_K         = 11    # minimizer k-mer (BLAST word_size equivalent)
MM2_W         = 5     # minimizer window (small → sensitive on tiny ~1KB ref)
MM2_A         = 2     # match score
MM2_B         = 3     # mismatch penalty (default=4; reduced for diverged tag tolerance)
MM2_GAP_OPEN  = "2,2" # two-piece gap open
MM2_GAP_EXT   = "1,1" # two-piece gap extension
MM2_MIN_SCORE = 15    # minimum alignment score

_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def _cigar_metrics(cigar: str) -> tuple[int, int, int]:
    """
    Parse CIGAR string and return (leading_clips, aln_query_len, aln_ref_len).

    leading_clips  — soft/hard clips before the alignment = 0-based query start
    aln_query_len  — bases consumed in the query by the alignment (M/I/=/X)
    aln_ref_len    — bases consumed in the reference by the alignment (M/D/N/=/X)
    """
    leading = 0
    qlen = 0
    rlen = 0
    seen = False
    for n_str, op in _CIGAR_RE.findall(cigar):
        length = int(n_str)
        if op in ("S", "H") and not seen:
            leading += length
        elif op in ("M", "=", "X"):
            seen = True
            qlen += length
            rlen += length
        elif op == "I":
            seen = True
            qlen += length
        elif op in ("D", "N"):
            seen = True
            rlen += length
    return leading, qlen, rlen


def run_minimap2_search(
    tags_fasta: Path,
    reads_inputs: list[Path],
    out_sam: Path,
    threads: int = 4,
) -> None:
    """
    Map WGS reads to tag library (tiny ~1KB reference) with minimap2.

    All reads files are passed as separate inputs (single-end mode — we don't
    need paired-end insert size, only which reads contain a tag).
    --secondary=yes -N 100: allow multiple hits per read (one per tag).

    Output is piped through `samtools view -F 4 -h` to keep only mapped reads
    plus the SAM header. This prevents a 30M-line file of unmapped reads and
    reduces parse time from ~20s to <1s.
    """
    mm2_cmd = [
        "minimap2",
        "-a",                   # SAM output
        "-k", str(MM2_K),
        "-w", str(MM2_W),
        "-A", str(MM2_A),
        "-B", str(MM2_B),
        "-O", MM2_GAP_OPEN,
        "-E", MM2_GAP_EXT,
        "-s", str(MM2_MIN_SCORE),
        "--secondary=yes",
        "-N", "100",            # allow up to 100 secondary alignments per read
        "-t", str(threads),
        str(tags_fasta),
    ] + [str(r) for r in reads_inputs]

    filter_cmd = ["samtools", "view", "-F", "4", "-h"]

    with open(out_sam, "w") as f:
        mm2 = subprocess.Popen(
            mm2_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL
        )
        filt = subprocess.Popen(
            filter_cmd, stdin=mm2.stdout, stdout=f, stderr=subprocess.DEVNULL
        )
        mm2.stdout.close()
        filt.wait()
        mm2.wait()
        if mm2.returncode != 0:
            raise subprocess.CalledProcessError(mm2.returncode, mm2_cmd)
        if filt.returncode != 0:
            raise subprocess.CalledProcessError(filt.returncode, filter_cmd)


def parse_mm2_hits_and_extract_tracts(
    sam_path: Path,
    tag_info: dict[str, str],
    hits_tsv: Path,
    tracts_fasta: Path,
    min_aln_len: int = 5,
) -> tuple[int, int]:
    """
    Single-pass SAM parsing: apply anchor filter, extract 20-nt tracts inline.

    Reference = tag library, query = WGS read.

    Anchor filter (on reference = tag coordinates):
      START tag: alignment must cover the 5' anchor → POS <= ANCHOR
      END tag:   alignment must reach the 3' anchor → ref_end >= TAG_LEN - ANCHOR + 1

    Tract extraction uses the SAM SEQ field directly.
    For minus-strand alignments (FLAG & 16), SAM already stores the reverse
    complement of the original read, so query coordinates work without flipping.

    Writes hits.tsv and tracts.fa in identical format to the BLAST path.
    Returns (n_hits, n_tracts).
    """
    hits = []
    tract_records = []

    with open(sam_path) as f:
        for line in f:
            if line.startswith("@"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 10:
                continue

            flag  = int(parts[1])
            rname = parts[2]       # reference = tag name
            pos   = int(parts[3])  # 1-based alignment start on tag
            cigar = parts[5]
            seq   = parts[9]       # query seq (RC'd by minimap2 if FLAG&16)

            # Skip unmapped and supplementary
            if flag & 4 or flag & 2048:
                continue
            if rname == "*" or cigar == "*" or seq == "*":
                continue

            tag_type = tag_info.get(rname, "")
            if not tag_type:
                continue

            leading, aln_qlen, aln_rlen = _cigar_metrics(cigar)

            if aln_rlen < min_aln_len:
                continue

            # 1-based inclusive end on the tag (reference)
            ref_end = pos + aln_rlen - 1

            # Anchor filter — same logic as parse_blast_hits but on ref coords
            if tag_type == "START" and pos > ANCHOR:
                continue
            if tag_type == "END" and ref_end < TAG_LEN - ANCHOR + 1:
                continue

            # Query (read) coordinates for tract extraction
            r_start = leading            # 0-based start in SEQ
            r_end   = leading + aln_qlen # 0-based exclusive end in SEQ

            tract = extract_tract_from_read(seq, r_start, r_end, tag_type)
            if tract is None:
                continue

            hit_idx = len(hits)
            hits.append({
                "read_id":  parts[0],
                "tag_name": rname,
                "tag_type": tag_type,
                "r_start":  r_start,
                "r_end":    r_end,
                "strand":   "-" if (flag & 16) else "+",
            })
            tract_records.append((f"tract_{hit_idx}|{parts[0]}|{rname}", tract))

    write_hits_tsv(hits, hits_tsv)

    with open(tracts_fasta, "w") as f:
        for name, seq in tract_records:
            f.write(f">{name}\n{seq}\n")

    return len(hits), len(tract_records)
