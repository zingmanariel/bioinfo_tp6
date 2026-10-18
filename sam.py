"""SAM/BAM reading: header, alignments, CIGAR parsing."""

import re


# A SAM file is tab-separated text. Lines starting with '@' are the header
# (skipped); every other line is one read, with these columns (1-based):
#
#   1 QNAME  2 FLAG  3 RNAME  4 POS  5 MAPQ  6 CIGAR  7 RNEXT  8 PNEXT
#   9 TLEN  10 SEQ  11 QUAL  12+ optional tags
#
# Example (SEQ and QUAL shortened):
#
#   @HD  VN:1.6  SO:coordinate
#   @SQ  SN:chr1  LN:248956422
#   K00171:229:HCCKGBBXX:8:1209:24170:28252  99  chr1  209625384  60  101M  =  209625461  178  TTGGGCTCCT...  FENLLLLLJM...  MC:Z:101M ...
#
# We keep POS (col 4, int), CIGAR (col 6) and SEQ (col 10) -> indexes 3, 5, 9
# after split('\t'), since Python counts from 0.
FLAG_COL = 1
CHROM_COL = 2
POS_COL = 3
CIGAR_COL = 5
SEQ_COL = 9


def read_sam(path):
    """Alignments of a SAM file as a list of dicts (header lines are skipped)."""
    reads = []
    with open(path) as f:
        for line in f:
            if line.startswith('@'):
                continue
            cols = line.rstrip('\n').split('\t')
            reads.append({
                'flag': int(cols[FLAG_COL]),
                'chrom': cols[CHROM_COL],
                'pos': int(cols[POS_COL]),     # 1-based, first covered position
                'cigar': cols[CIGAR_COL],
                'seq': cols[SEQ_COL],
            })
    return reads


# CIGAR operations that move along the reference (the others, S and I, only
# consume bases of the read).
REFERENCE_OPS = ('M', 'D', 'N', '=', 'X')


def parse_cigar(cigar):
    """'96M5S' -> [(96, 'M'), (5, 'S')]"""
    return [(int(n), op) for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar)]


def covered_length(cigar):
    """How many reference positions a CIGAR covers (M, D, N, =, X count; S, I don't)."""
    return sum(n for n, op in parse_cigar(cigar) if op in REFERENCE_OPS)


DUPLICATE_FLAG = 1024   # bit set by MarkDuplicates


def is_duplicate(read):
    return bool(read['flag'] & DUPLICATE_FLAG)
