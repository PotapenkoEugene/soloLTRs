"""
tags.py — Extract START/END tags and 20-nt tracts from LTR-RT representative paralogs.

For each full-length paralog (Cossu Dataset S1 or EDTA-derived), extract:
  - START tag: first TAG_LEN nt of the paralog (= beginning of 5' LTR)
  - END tag  : last  TAG_LEN nt of the paralog (= end of 3' LTR)

These tags are used to search WGS reads (RepeatMasker/BLAST step).
When a read matches a tag, a 20-nt TRACT is extracted from the read:
  - For START tag match at position p in read:
      tract = read[p - FLANK : p + ANCHOR]   (15 nt upstream + 5 nt from tag start)
  - For END tag match ending at position q in read:
      tract = read[q - ANCHOR : q + FLANK]   (5 nt from tag end + 15 nt downstream)

Tracts are then mapped to the paralog library (BWA ALN):
  M (mapped)   = tract comes from internal region of a complete element
  U (unmapped) = tract comes from genomic flanking (solo LTR or element boundary)

Formula: S/C = U/M - 1   (Cossu et al. 2017, GBE)
"""

from pathlib import Path

TAG_LEN = 50     # length of START/END tag
ANCHOR  = 5      # nt taken from within the tag for the tract
FLANK   = 15     # nt taken from the read flanking beyond the tag
TRACT_LEN = ANCHOR + FLANK   # = 20 nt


def load_fasta(path: Path) -> list[tuple[str, str]]:
    """Load FASTA → list of (name, seq) preserving order."""
    records: list[tuple[str, str]] = []
    name = None
    parts: list[str] = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(parts).upper()))
                name = line[1:].split()[0]
                parts = []
            elif name is not None:
                parts.append(line)
    if name is not None:
        records.append((name, "".join(parts).upper()))
    return records


def write_fasta(records: list[tuple[str, str]], path: Path, wrap: int = 60) -> None:
    with open(path, "w") as f:
        for name, seq in records:
            f.write(f">{name}\n")
            for i in range(0, len(seq), wrap):
                f.write(seq[i : i + wrap] + "\n")


def extract_tags(paralogs: list[tuple[str, str]]) -> list[tuple[str, str, str]]:
    """
    For each paralog, return (tag_name, tag_type, tag_seq).
    tag_type: "START" or "END"
    tag_name: e.g. "Arabidopsis_1_START" or "Arabidopsis_1_END"
    """
    tags = []
    for name, seq in paralogs:
        if len(seq) < TAG_LEN:
            continue
        tags.append((f"{name}_START", "START", seq[:TAG_LEN]))
        tags.append((f"{name}_END",   "END",   seq[-TAG_LEN:]))
    return tags


def build_tag_fasta(paralogs_path: Path, out_tags_path: Path) -> list[tuple[str, str, str]]:
    """
    Read paralogs FASTA, extract START/END tags, write tag library FASTA.
    Returns list of (tag_name, tag_type, tag_seq).
    """
    paralogs = load_fasta(paralogs_path)
    tags = extract_tags(paralogs)
    records = [(tag_name, tag_seq) for tag_name, _, tag_seq in tags]
    write_fasta(records, out_tags_path)
    return tags


def extract_tract_from_read(
    read_seq: str,
    match_start: int,
    match_end: int,
    tag_type: str,
) -> str | None:
    """
    Given a read and a tag match position, extract the 20-nt tract.

    For START tag: tract = read[match_start - FLANK : match_start + ANCHOR]
                         = 15 nt before the tag start + first 5 nt of the tag
    For END tag:   tract = read[match_end - ANCHOR : match_end + FLANK]
                         = last 5 nt of the tag + 15 nt after the tag end

    Returns None if there are not enough flanking nucleotides.
    """
    if tag_type == "START":
        lo = match_start - FLANK
        hi = match_start + ANCHOR
    elif tag_type == "END":
        lo = match_end - ANCHOR
        hi = match_end + FLANK
    else:
        raise ValueError(f"Unknown tag_type: {tag_type!r}")

    if lo < 0 or hi > len(read_seq):
        return None
    return read_seq[lo:hi]
