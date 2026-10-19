"""Ejercicio 61b) cargue el BAM en el IGV y trate de obtener con el mismo las métricas obtenidas
en 61a, compare los resultados.

Region de interes: chr1:209,624,476-209,628,968
"""


import os

from depth import THRESHOLDS, coverage, mean_sd
from paths import FASTA, FIGURES, SAM
from pileup import Pileup
from plots import plot_depth
from pileup import read_sam
from variants import alts_at

def run():
    pileup = Pileup.from_files(SAM, FASTA)
    print(alts_at(pileup, 209625908))
    print(alts_at(pileup, 209625710))
