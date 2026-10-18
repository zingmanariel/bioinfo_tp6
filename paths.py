"""Paths to the data files and output folders used by the TP."""

import os

ROOT = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(ROOT, 'data')
FIGURES = os.path.join(ROOT, 'figures')

SAM = os.path.join(DATA, 'modulo6.sam')
FASTA = os.path.join(DATA, 'Homo_sapiens_assembly38.chr1.fasta.gz')
