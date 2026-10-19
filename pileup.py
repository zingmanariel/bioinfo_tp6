"""Loading of the data (SAM reads, reference fasta) and the Pileup that holds both.

Pileup: the reference and the reads loaded together, one counter per position.

    chrom   chromosome name (e.g. chr1)
    start   first position covered by any read (1-based, like POS in the SAM)
    ref     reference bases for the window [start, start + len(ref) - 1]
    pile    pile[i] = {base: count} for position start + i; '*' counts deletions
"""

import gzip
import re
from collections import Counter
from itertools import islice


# ---------------------------------------------------------------------------
# SAM
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Reference fasta
# ---------------------------------------------------------------------------

LINE = 60   # bases per line in the fasta (see the .fai)


def fetch_reference(path, start, end):
    """Reference bases for the 1-based inclusive range [start, end], upper-cased.

    The file is a single chromosome, so the base at 0-based offset k sits in
    line 1 + k // LINE (line 0 is the '>chr1' header). Only those lines are kept.
    """
    k0, k1 = start - 1, end            # 0-based, end exclusive
    first, last = 1 + k0 // LINE, 1 + (k1 - 1) // LINE
    with gzip.open(path, 'rt') as f:
        lines = islice(f, first, last + 1)
        seq = ''.join(line.strip() for line in lines).upper()
    offset = k0 - (first - 1) * LINE
    return seq[offset:offset + (end - start + 1)]


# ---------------------------------------------------------------------------
# Pileup
# ---------------------------------------------------------------------------
REF_OPS = 'M=X'     # consume read and reference, base present
DEL = '*'


class Pileup:
    def __init__(self, chrom, start, ref, pile):
        self.chrom = chrom
        self.start = start
        self.ref = ref
        self.pile = pile

    @classmethod
    def from_files(cls, sam_path, fasta_path, skip_duplicates=True):
        """The window always spans all the reads; with skip_duplicates the reads
        marked as duplicates (flag 1024) are not counted, as IGV and GATK do."""
        # load reads
        reads = read_sam(sam_path)

        # calculate the window covered by the reads
        start = min(read['pos'] for read in reads)
        end = 0
        for read in reads:
            read_end = read['pos'] + covered_length(read['cigar']) - 1
            end = max(end, read_end)
        window_size = end - start + 1

        # one empty counter per position of the window where we can add each base separately
        pile = [Counter() for _ in range(window_size)]
        for read in reads:
            if skip_duplicates and is_duplicate(read):
                continue
            cls._add_read(pile, read, start)

        # load the reference in the covered window
        ref = fetch_reference(fasta_path, start, end)
        return cls(reads[0]['chrom'], start, ref, pile)

    @staticmethod
    def _add_read(pile, read, start):
        """Add one read to the counters, walking its CIGAR."""
        seq = read['seq']
        window_index = read['pos'] - start    # where the read is in the window
        seq_index = 0                         # where we are in the read's bases

        for length, op in parse_cigar(read['cigar']):
            if op in REF_OPS:
                for k in range(length):
                    base = seq[seq_index + k]
                    counter = pile[window_index + k]
                    counter[base] += 1
                window_index += length
                seq_index += length
            elif op == 'D':
                for k in range(length):
                    counter = pile[window_index + k]
                    counter[DEL] += 1
                window_index += length
            elif op == 'N':
                window_index += length
            elif op in 'IS':
                seq_index += length           # bases in the read, not in the reference
            # H and P consume nothing

    def __len__(self):
        return len(self.pile)

    def counts_at(self, pos):
        return self.pile[pos - self.start]

    def ref_at(self, pos):
        return self.ref[pos - self.start]

    def depth_at(self, pos):
        return sum(self.counts_at(pos).values())

    def depths(self):
        return [sum(c.values()) for c in self.pile]
