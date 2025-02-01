import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt


def proportions_for_barplot(adata, ct_key="cell type"):
    """Create df with proportions of CNV clusters within each cell type.

    Args:
        adata: AnnData object
        ct_key: name of the cell type column

    Return:
        props: df with proportions
    """

    sizes = adata.obs.groupby([ct_key, "cnv_leiden"]).size()
    props = pd.DataFrame(sizes.groupby(ct_key).apply(lambda x: 100 * x/ x.sum()).reset_index(level=0, drop=True).reset_index())
    props = props.pivot(columns=ct_key, index="cnv_leiden").T
    props.index = props.index.droplevel(0)

    return props

def proportions_gt_for_barplot(adata, tp53_key):
    """Create df with proportions of TP53 mutations within each CNV cluster.

    Args:
        adata: AnnData object
        tp53_key: column name for TP53 mutation information

    Returns:
        props: df with proportions
    """

    sizes = adata.obs.groupby(["cnv_leiden", tp53_key]).size()
    props = pd.DataFrame(sizes.groupby("cnv_leiden").apply(lambda x: 100 * x/ x.sum()).reset_index(level=0, drop=True).reset_index())
    props = props.pivot(columns="cnv_leiden", index=tp53_key).T
    props.index = props.index.droplevel(0)

    return props

def plot_cluster_proportions(cluster_props, 
                             cluster_palette=None,
                             xlabel_rotation=0,
                            figsize=(8,8),
                            filename=None): 
    """Plot barplots with the cluster proportion df.

    Args:
       cluster_props: df with proportions
    """

    fig, ax = plt.subplots(dpi=300, figsize=figsize)
    fig.patch.set_facecolor("white")

    cmap = None
    if cluster_palette is not None:
        cmap = sns.palettes.blend_palette(
            cluster_palette, 
            n_colors=len(cluster_palette), 
            as_cmap=True)

    cluster_props.plot(
        kind="bar", 
        stacked=True, 
        ax=ax, 
        legend=None, 
        colormap=cmap
    )

    #reorder legend labels
    handles, labels = ax.get_legend_handles_labels()
    labels = list(reversed(labels))
    handles = list(reversed(handles))
    ax.legend(handles, labels, bbox_to_anchor=(1.01, 1), frameon=False, title="Cluster")
    sns.despine(fig, ax)
    ax.tick_params(axis="x", rotation=xlabel_rotation)
    #ax.set_xlabel(cluster_props.index.name.capitalize())
    ax.set_xlabel("")
    ax.set_ylabel("Proportion")
    ax.grid(False)
    fig.tight_layout()
    fig.savefig(filename, dpi=300)

    return fig

