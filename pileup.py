"""Pileup: the reference and the reads loaded together, one counter per position.

    start   first position covered by any read (1-based, like POS in the SAM)
    ref     reference bases for the window [start, start + len(ref) - 1]
    pile    pile[i] = {base: count} for position start + i; '*' counts deletions
"""

from collections import Counter

import fasta
from sam import covered_length, is_duplicate, parse_cigar, read_sam

REF_OPS = 'M=X'     # consume read and reference, base present
DEL = '*'


class Pileup:
    def __init__(self, start, ref, pile):
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
        ref = fasta.fetch(fasta_path, start, end)
        return cls(start, ref, pile)

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
