"""Shared plotting utilities used across figure-generating scripts."""

import os

import matplotlib.pyplot as plt

FIGURE_DPI = 300


def significance_stars(p):
    """Convert a p-value to a standard star notation for plots and logs."""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def save_figure(script, path, bbox_inches=None):
    """Save the current matplotlib figure, close it, and log the output path."""
    kwargs = {"dpi": FIGURE_DPI}
    if bbox_inches is not None:
        kwargs["bbox_inches"] = bbox_inches
    plt.savefig(path, **kwargs)
    plt.close()
    print(f"[{script}] Saved: {os.path.abspath(path)}")
