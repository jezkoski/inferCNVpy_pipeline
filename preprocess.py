import pandas as pd
import numpy as np
import scanpy as sc


def read_data(path: str):
    """Read anndata object from the path.
    
    Args:
        path (str): filepath

    Returns:
        adata: AnnData object

    """

    adata = sc.read_h5ad(path)

    return adata



def annot_biomart(adata, biomart_file_path):
    """Annotate adata object with biomart gene position. 

    Adds columns from biomart file to the adata object.

    Args:
        adata: AnnData object
        biomart_file_path: path to biomart file that has columns 
                           ["chromosome_name", "ensembl_gene_id", "start_position", "end_position"]

    """

    biomartannot = pd.read_csv(biomart_file_path, sep=" ")
    #add 'chr' to the chromosome_name because inferCNVpy doesn't run otherwise
    biomartannot["chromosome_name"] = "chr" + biomartannot["chromosome_name"]
    biomartannot.index = biomartannot["ensembl_gene_id"]

    #index of adata and file has to be the same
    adata.var["gene_name"] =adata.var_names
    adata.var_names = adata.var["gene_ids"]

    adata.var["chromosome"] = biomartannot["chromosome_name"]
    adata.var["start"] = biomartannot["start_position"]
    adata.var["end"] = biomartannot["end_position"]
    adata.var_names = adata.var["gene_name"]



def clean_adata(adata, chromosomes=None):
    """Clean up chromosome and position columns for the cnv analysis.

    Args:
        adata: AnnData object
        chromosomes: list of chromosomes to remove

    Returns:
        adata: cleaned AnnData object with problematic chromosomes removed
    """

    adata = adata[:,~adata.var["chromosome"].isna()]
    adata.var["ensg"] = adata.var["gene_ids"] #rename this one too (same as infercnvpy docs)

    if chromosomes is None:
        #if no chromosome file provided, remove these chromosomes from the file

        exclude_chromosomes = [ 'GL000009.2','GL000194.1','GL000195.1',
                                'GL000218.1', 'GL000219.1', 'KI270711.1', 'KI270721.1',
                                'KI270726.1', 'KI270727.1', 'KI270728.1', 'KI270731.1',
                                'KI270734.1', 'MT','chrGL000009.2','chrGL000194.1',
                                'chrGL000195.1', 'chrGL000218.1', 'chrGL000219.1', 'chrKI270711.1',
                                'chrKI270721.1', 'chrKI270726.1', 'chrKI270727.1', 'chrKI270728.1', 
                                'chrKI270731.1','chrKI270734.1', 'chrMT']

        adata = adata[:,~adata.var["chromosome"].isin(exclude_chromosomes)]
    else:
        adata = adata[:, ~adata.var["chromosome"].isin(chromosomes)]

    return adata


def remove_samples(adata, samples: list):
    """Remove samples from the adata object

    Args:
        adata: AnnData object
        samples: list of sample names to be removed
    Return:
        adata_subset: subseted adata object
    """

    adata_subset = adata.copy()
    adata_subset = adata_subset[~adata_subset.obs["sample"].isin(samples)]

    return adata_subset


