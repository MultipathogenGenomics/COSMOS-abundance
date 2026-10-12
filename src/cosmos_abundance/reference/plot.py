# Plot reference-set Mash trees.

import matplotlib
matplotlib.use("Agg")

from matplotlib import pyplot as plt
import numpy as np
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform

def plot_mash_tree(ids: list[str], dist_matrix: np.ndarray, representative_ids: list[str], mash_distance_threshold: float, output_file: str, labels_by_id: dict[str, str]) -> None:
    if len(ids) < 2:
        fig, ax = plt.subplots(figsize=(6, 2))
        if len(ids) == 1:
            ax.text(0.02, 0.5, labels_by_id[ids[0]], ha="left", va="center", fontweight="bold")
        else:
            ax.text(0.02, 0.5, "No sequences to plot", ha="left", va="center")
        ax.set_axis_off()
        fig.tight_layout()
        fig.savefig(output_file, dpi=300)
        plt.close(fig)
        return
    condensed_distance = squareform(dist_matrix, checks=False)
    linkage_matrix = sch.linkage(condensed_distance, method="complete", optimal_ordering=True)
    fig, ax = plt.subplots(figsize=(8, max(6, len(ids) * 0.16)))
    labels = [labels_by_id[seq_id] for seq_id in ids]
    dendrogram = sch.dendrogram(
        linkage_matrix,
        labels=labels,
        ax=ax,
        leaf_rotation=0,
        leaf_font_size=8,
        orientation="left",
        above_threshold_color="black",
        color_threshold=0,
    )
    representative_ids = set(representative_ids)
    for label, leaf_index in zip(ax.get_yticklabels(), dendrogram["leaves"]):
        if ids[leaf_index] in representative_ids:
            label.set_fontweight("bold")
    ax.axvline(mash_distance_threshold, color="black", linestyle="--", linewidth=1, label=f"threshold: {mash_distance_threshold:g}")
    ax.set_xlabel("Mash-like distance (complete linkage)")
    ax.legend(frameon=False, loc="upper left", fontsize=8)
    ax.grid(axis="x", linestyle="--", linewidth=0.6, alpha=0.35)
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.tight_layout()
    fig.savefig(output_file, dpi=300)
    plt.close(fig)
