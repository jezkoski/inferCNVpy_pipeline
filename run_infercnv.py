import infercnvpy as cnv
from anndata import AnnData
import scanpy as sc


def run_infercnv(adata: AnnData, reference_key=None, reference_cat=None):
	"""Run inferCNVpy and provide normal cell types

	Args:
		adata: AnnData object
		reference_key: column name with reference cell type names
		reference_cat: cell type names which are used as reference
	"""

	if reference_key is None:

		cnv.tl.infercnv(
			adata,
			window_size=250
		)
	else:
		cnv.tl.infercnv(
			adata,
			reference_key=reference_key,
			reference_cat=reference_cat,
			window_size=250
		)


def cluster_by_cnvprofile(adata: AnnData):
	"""Cluster the cells by CNV profiles with leiden clustering.

	Args:
        	adata: AnnData object
	"""

	cnv.tl.pca(adata)
	cnv.pp.neighbors(adata)
	cnv.tl.leiden(adata)


def umap_and_score(adatA: AnnData):
	"""Create an umap and score the cnvs.
	"""

	cnv.tl.umap(adata)
	cnv.tl.cnv_score(adata)

