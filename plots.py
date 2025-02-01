import pandas as pd
import numpy as np
import matplotpĺib.pyplot as plt
from matplotlib.pyplot import rc_context
import scanpy as sc
import infercnvpy as cnv


sc.settings.set_figure_params(figsize=(5,5), fontsize=20)

def plot_chromosome_heatmap(adata, groupby="cell_type",figsize=(15,23), dendrogram=False, filename=None):
    """Plot heatmap by chromosome (columns) and groups on rows.

    Args:
        adata: AnnData object
        groupby: grouping for rows in heatmap, eg. cell_type or cnv_cluster
        filename: filename to give to the saved plot
        dendrogram: show dendrogram, default False
        figsize: figure size as a tuple (10,10)
    """

    cnv.pl.chromosome_heatmap(adata, groupby=groupby, dendrogram=dendrogram, save=filename, figsize=figsize)


def plot_chromosome_heatmap_summary(adata, groupby="cnv_leiden",figsize=(10,15), filename=None):
    """Plot heatmap summary by chromosome (columns) and groups on rows. Shows CNV segmentation like plot.

    Args:
        adata: AnnData object
        groupby: grouping for rows in heatmap, eg. cell_type or cnv_cluster
        figsize: figure size
        filename: filename to give to the saved plot
    """

    cnv.pl.chromosome_heatmap_summary(adata, groupby=groupby, figsize=figsize, save=filename)


def plot_3_cnv_umaps(adata, filename=None, colors=["cnv_leiden", "cnv_score", "cell_type"]):
    """Plots 3 UMAPs with infercnvpy tool.

    Args:
        adata: AnnData object
        filename: name for saved figure
        colors: list of colors to use in the 3 UMAPS eg. ["cnv_leiden", "cnv_score", "cell_type"]

    """

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11, 11))
    ax4.axis("off")
    cnv.pl.umap(
        adata,
        color=colors[0],
        legend_loc="on data",
        legend_fontoutline=2,
        ax=ax1,
        show=False,
    )
    cnv.pl.umap(adata, color=colors[1], ax=ax2, show=False)
    cnv.pl.umap(adata, color=colors[2], ax=ax3)
    fig.savefig(filename, dpi=300)


def plot_3_scanpy_umaps(adata, filename=None, colors=["cnv_leiden", "cnv_score", "cell_type"]):
    """Plots 3 UMAPs with scanpy.

    Args:
        adata: AnnData object
        filename: name for saved figure
        colors: list of colors to use in the 3 UMAPS eg. ["cnv_leiden", "cnv_score", "cell_type"]

    """

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 11), gridspec_kw={"wspace": 0.5})
    ax4.axis("off")
    sc.pl.umap(adata, color=colors[0], ax=ax1, show=False)
    sc.pl.umap(adata, color=colors[1], ax=ax2, show=False)
    sc.pl.umap(adata, color=colors[2], ax=ax3)

    fig.savefig(filename, dpi=300)

def plot_single_umap(adata, color, filename=None):
    """Plot a single UMAP with scanpy

    Args:
        adata: AnnData object
        color: color to use in the UMAP
        filename: name of the figure

    """

    sc.pl.umap(adata, color=color, size=40, save=filename)


