"""Reference genome loading (Homo_sapiens_assembly38.chr1.fasta.gz)."""

import gzip
from itertools import islice

LINE = 60   # bases per line in the fasta (see the .fai)


def fetch(path, start, end):
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
