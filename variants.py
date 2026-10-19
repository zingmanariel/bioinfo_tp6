"""Candidate variant detection, depth-based filtering and VCF writing (6.2b)."""

# Filter thresholds (6.2b-ii)
MIN_DEPTH = 10           # minimum reads covering the position
MIN_ALT_HET = 5          # minimum reads supporting a heterozygous alt
MIN_ALT_HOM = 5          # minimum reads supporting a homozygous alt
HET_MIN_FRACTION = 0.20  # alt fraction range that counts as heterozygous: [0.20, 0.80)
HOM_MIN_FRACTION = 0.80  # alt fraction from which it counts as homozygous


def alts_at(pileup, pos):
    """Return the bases that are not the reference at a given position."""
    ref = pileup.ref_at(pos)
    counts = pileup.counts_at(pos)
    return {base: count for (base, count) in counts.items() if base != ref}


def candidate_variants(pileup):
    """One dict per (position, alt) found in the window.

    pos       1-based position
    ref       reference base
    alt       observed base that differs from the reference ('*' = deletion)
    depth     reads covering the position
    alt_depth reads supporting the alt
    """
    variants = []
    for pos in range(pileup.start, pileup.start + len(pileup)):
        depth = pileup.depth_at(pos)
        ref = pileup.ref_at(pos)
        for alt, alt_depth in alts_at(pileup, pos).items():
            variants.append({
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'depth': depth,
                'alt_depth': alt_depth,
            })
    return variants


def genotype(depth, alt_depth):
    """'1/1' (homozygous alt), '0/1' (heterozygous) or None if discarded."""
    if depth < MIN_DEPTH:
        return None

    alt_fraction = alt_depth / depth
    if alt_fraction >= HOM_MIN_FRACTION and alt_depth >= MIN_ALT_HOM:
        return '1/1'
    if HET_MIN_FRACTION <= alt_fraction < HOM_MIN_FRACTION and alt_depth >= MIN_ALT_HET:
        return '0/1'
    return None


def write_vcf(variants, chrom, path):
    """Write the variants as a minimal VCF. Each variant needs pos, ref, alt,
    depth, alt_depth and genotype (the dicts from variants.candidate_variants)."""
    with open(path, 'w') as f:
        f.write('##fileformat=VCFv4.2\n')
        f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for v in variants:
            fraction = v['alt_depth'] / v['depth']
            info = f"DP={v['depth']};AD={v['alt_depth']};AF={fraction:.3f};GT={v['genotype']}"
            f.write(f"{chrom}\t{v['pos']}\t.\t{v['ref']}\t{v['alt']}\t.\tPASS\t{info}\n")
