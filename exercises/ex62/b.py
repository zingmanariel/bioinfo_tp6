"""Ejercicio 62b): En base a los archivos SAM y FASTA y lo analizado en el ejercicio 61a, realice
un programa que determine las variantes y las escriba en un archivo simil VCF. Para ello
programe los siguientes módulos y luego combinelos
i) realice un programa que determine todas las potenciales variantes y para cada una guarde:
a) la profundidad de lecturas en “esa” posición del genoma de referencia.
b) la profundidad de la potencial variante observada

ii) realice un código que dados como parámetros umbrales de Profundidad de la posición y
profundidad de la variante (para variantes homo y heterocigotas) filtre las diferencias y escriba
el VCF resultante con las variantes.
"""

from paths import FASTA, SAM, VCF_OUT
from pileup import Pileup
from variants import candidate_variants, genotype, write_vcf


def i():
    variants = candidate_variants(Pileup.from_files(SAM, FASTA))
    print(f'{len(variants)} candidate variants')
    for v in variants:
        print(f"{v['pos']}  {v['ref']}>{v['alt']}  depth {v['depth']}  alt {v['alt_depth']}")


def ii():
    pileup = Pileup.from_files(SAM, FASTA)
    kept = []
    for v in candidate_variants(pileup):
        v['genotype'] = genotype(v['depth'], v['alt_depth'])
        if v['genotype']:
            kept.append(v)
    write_vcf(kept, pileup.chrom, VCF_OUT)
    print(f'{len(kept)} variants written to {VCF_OUT}')
    for v in kept:
        print(f"{v['pos']}  {v['ref']}>{v['alt']}  {v['genotype']}  depth {v['depth']}  alt {v['alt_depth']}")
