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


def add_significance_legend(ax, loc="upper right"):
    """Add a corner legend explaining Mann-Whitney significance stars."""
    legend_text = (
        "Mann-Whitney U test\n"
        "* p < 0.05\n"
        "** p < 0.01\n"
        "*** p < 0.001\n"
        "ns not significant"
    )
    coordinates = {
        "upper right": (0.98, 0.98, "right", "top"),
        "upper left": (0.02, 0.98, "left", "top"),
    }
    x, y, ha, va = coordinates.get(loc, coordinates["upper right"])
    ax.text(
        x, y, legend_text,
        transform=ax.transAxes,
        ha=ha, va=va, fontsize=8,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85, edgecolor="0.8"),
    )


def save_figure(script, path, bbox_inches=None):
    """Save the current matplotlib figure, close it, and log the output path."""
    kwargs = {"dpi": FIGURE_DPI}
    if bbox_inches is not None:
        kwargs["bbox_inches"] = bbox_inches
    plt.savefig(path, **kwargs)
    plt.close()
    print(f"[{script}] Saved: {os.path.abspath(path)}")
