"""Ejercicio 61a) Obtenga a partir del archivo SAM las siguientes métricas
i) Profundidad para cada una de las posiciones del genoma de referencia (realice el gŕafico
correspondiente)
ii) Profundidad promedio y desvío estándar
iii) Promedio y SD del largo de las lecturas mapeadas
iv) Cobertura (en %) del genoma target para un dado umbral de profundidad.

The "target" is the window spanned by the reads (first to last covered
position); positions inside it that no read covers count as depth 0.
"""

import os

from depth import THRESHOLDS, coverage, mean_sd
from paths import FASTA, FIGURES, SAM
from pileup import Pileup
from plots import plot_depth
from pileup import read_sam


def i():
    pileup = Pileup.from_files(SAM, FASTA)
    depths = pileup.depths()
    end = pileup.start + len(pileup) - 1
    print(f'window: chr1:{pileup.start}-{end} ({len(pileup)} positions)')
    print(f'positions with depth 0: {depths.count(0)}')
    path = os.path.join(FIGURES, '61a_depth.png')
    plot_depth(pileup.start, depths, path)
    print(f'figure: {path}')


def ii():
    depths = Pileup.from_files(SAM, FASTA).depths()
    mean, sd = mean_sd(depths)
    print(f'depth over the window ({len(depths)} positions): mean {mean:.2f}, SD {sd:.2f}')


def iii():
    lengths = [len(r['seq']) for r in read_sam(SAM)]
    mean, sd = mean_sd(lengths)
    print(f'{len(lengths)} reads: length mean {mean:.2f}, SD {sd:.2f} '
          f'(min {min(lengths)}, max {max(lengths)})')


def iv():
    depths = Pileup.from_files(SAM, FASTA).depths()
    for t in THRESHOLDS:
        print(f'coverage at depth >= {t:>3}: {coverage(depths, t):6.2f}%')
