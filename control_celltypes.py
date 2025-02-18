import scanpy as sc
import pandas as pd
from typing import Dict


def set_control_celltypes(adata, healthy_samples, ct_key_new):
    """Add (non-malignant) tag to the cell types for control samples.

    Args:
        adata: AnnData object
        healthy_samples: list of the names of healthy samples
        ct_key_new: name of the cell type column used
    """

    temp = pd.DataFrame(adata.obs[["sample",ct_key_new]])

    ct_annotation = {'B cells': 'B cells (non-malignant)' ,
     'DC cells': 'DC cells (non-malignant)',
     'EMP': 'EMP (non-malignant)',
     'Early eryth prog': 'Early eryth prog (non-malignant)',
     'Granulopoietic cells': 'Granulopoietic cells (non-malignant)',
     'HSCs & MPPs': 'HSCs & MPPs (non-malignant)',
     'LMP': 'LMP (non-malignant)',
     'Late eryth prog': 'Late eryth prog (non-malignant)',
     'Mesenchymal cells': 'Mesenchymal cells (non-malignant)',
     'MkP': 'MkP (non-malignant)',
     'Monocytopoietic cells': 'Monocytopoietic cells (non-malignant)',
     'NK cells': 'NK cells (non-malignant)',
     'Neutrophils': 'Neutrophils (non-malignant)',
     'Plasma cells': 'Plasma cells (non-malignant)',
     'T cells': 'T cells (non-malignant)'}

    temp["ct_cnv"] = pd.DataFrame(temp.loc[temp["sample"].isin(healthy_samples), ct_key_new].map(ct_annotation))
    temp["test"]= (temp["ct_cnv"]== np.nan).astype(str)
    temp["ct_cnv"] = temp["ct_cnv"].astype(str, errors="ignore")
    temp.loc[temp["ct_cnv"]=="nan", "ct_cnv"] = temp.loc[temp["ct_cnv"]=="nan", ct_key_new]

    adata.obs["ct_cnv"] = temp["ct_cnv"]


def add_samplename_to_celltypes(adata, healthy_samples, ct_key_new) -> Dict[str, str]:
    """Add sample name to cell type name and save it as separate column 'ct_ref'

    Args:
        adata: AnnData object
        healthy_samples:  list of healthy samples
        ct_key_new: name of the cell type column used

    Return:
        ct_annot: Dictionary with regular cell types as keys and control cell types as values
    """

    temp = pd.DataFrame(adata[adata.obs["ct_cnv"].isin(reference_cat)].obs[["sample","ct_cnv",ct_key_new]])
    temp = temp.loc[temp["sample"].isin(control_sample)]
    temp["ct_ref"] = temp["sample"].astype(str) + " " + temp["ct_cnv"].astype(str)
    temp["ct_ref"] = temp["ct_ref"].astype(str, errors="ignore")
    adata.obs["ct_ref"] = temp["ct_ref"]

    ct_annotation = dict(zip(temp[ct_key_new],temp["ct_ref"]))

    return ct_annot

