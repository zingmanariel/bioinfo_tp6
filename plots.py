"""Figures for the report."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plot_depth(start, depths, path):
    """Depth along the reference window; x axis is the chr1 position."""
    plt.figure(figsize=(12, 4))
    plt.plot(range(start, start + len(depths)), depths, linewidth=0.8)
    plt.ticklabel_format(axis='x', useOffset=False, style='plain')
    plt.xlabel('Position on chr1')
    plt.ylabel('Depth')
    plt.title('Depth per position')
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()
