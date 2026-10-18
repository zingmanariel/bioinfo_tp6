"""Per-position depth, mean/SD, read length stats and coverage (6.1a)."""

import numpy as np

THRESHOLDS = [1, 5, 10, 20, 30, 50, 100]


def mean_sd(values):
    """Mean and (population) standard deviation."""
    values = np.asarray(values, dtype=float)
    return values.mean(), values.std()


def coverage(depths, threshold):
    """% of positions whose depth is at least `threshold`.

    Example: depths = [0, 3, 12, 40], threshold = 10
             -> 2 of the 4 positions reach 10x -> 50.0
    """
    positions_total = len(depths)
    positions_covered = sum(1 for depth in depths if depth >= threshold)
    return 100 * positions_covered / positions_total
