import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

def get_ct_gene(adata, gene, X,  ct_key):
    """Run linear regression gene_expression ~ X per cell type.

    Args:
        adata: AnnData object
        gene: name of the gene whose expression against X to test
        X: column name for a numerical feature for which to perform linear regression eg. cvn_score
        ct_key: column name for the cell types
    """

    gene_adata = adata[:,adata.var_names == gene]
    results = pd.DataFrame(columns=["slope", "intercept", "r", "p", "std_err"])

    for ct in set(adata.obs[ct_key]):

        try:
            slope, intercept, r, p, std_err = gene_expression_cnvScore(gene_adata[gene_adata.obs[ct_key]==ct], gene, X, ct)
            results.loc[ct, [ "slope", "intercept", "r", "p", "std_err"]] = [slope, intercept, r, p, std_err ]
        except:
            results.loc[ct, [ "slope", "intercept", "r", "p", "std_err"]] = [np.nan, np.nan, np.nan, np.nan, np.nan]

    results.to_csv(os.path.join(path,f"{gene}_expr_{X}_lm.txt"), sep="\t")


def gene_expression_cnvScore(adata, gene, X, cell=None):
    """Performs linear regression gene_expression ~ cnv_score, and plots scatter plot.

    Args:
        adata: AnnData object
        gene: name of the gene whose expression against X to test
        X: column name for a numerical feature for which to perform linear regression eg. cnv_score

    Returns:
        (slope, intercept, r, p, std_err): return results from the linear regression
    """

    gene_expr = adata.X.toarray()
    X_values = adata.obs[X]
    try:
        slope, intercept, r, p, std_err = stats.linregress(gene_expr[:,0], X_values)
    except:
    	return np.nan
    plt.scatter(X_values, gene_expr[:,0])
    plt.ylabel(f"{gene} expression")
    plt.xlabel(f"{X}")
    plt.title(cell)
    if cell == None:
        plt.savefig(f"{gene}_expression_{X}_scatter.png", dpi=300)
    else:
        plt.savefig(f"{cell}_{gene}_expression_{X}_scatter.png",dpi=300)

    plt.close()

    return (slope, intercept, r, p,std_err)
