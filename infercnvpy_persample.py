from plots import *
from preprocess import *
from proportion_plots import *
from run_infercnvpy import *
import scanpy as sc
import pandas as pd
from anndata import AnnData
import argparse, sys


def read_and_process_adata(path: str) -> AnnData:
	"""Read the adata object and preprocess it for inferCNVpy analysis.

	!Assumes you have a file "biomart-locations.txt" in your folder

	Args:
		path: path to the h5ad file

	Returns:
		adata: cleaned and biomart-annotated AnnData object
	"""

	adata = read_data(path)
	annot_biomart(adata, "biomart-locations.txt")
	clean_adata(adata)

	return adata


def run_by_sample(adata: AnnData, samples: list, control: str, reference_key: str, reference_cat:str):
	"""Run inferCNV analysis for the samples chosen.

	Args:
		adata: AnnData object
		samples: list of samples the analysis is done to

	"""


	for sample in samples:
		sample_adata = adata[adata.obs["sample"].isin([sample, control])

		run_infercnv(sample_adata, reference_key, reference_cat)
		cluster_by_cnvprofile(sample_adata)
		umap_and_score(adata)
		# do visuals



def main():

	parser = argparse.ArgumentParser(description="inferCNVpy analysis for multiple samples.")

	parser.add_argument("-?", action="help", help=argparse.SUPPRESS_HELP, dest="help")
	parser.add_argument("-i", "--input", dest="input", action="store", help="Path to the input h5ad file.", metavar="data/adata.h5ad", required="True")
	parser.add_argument("-s", "--samples", nargs="+", dest="samples", action="store", help="List of samples the analysis is run for", metavar="-s sample1 sample2 sample3", required="True")

	args = parser.parse_args()

	if len(sys.argv)==1:
		parser.print_help()
		sys.exit()

	adata = read_and_process_adata(args.input)
	samples = args.samples

	run_by_sample(adata, samples)


if __name__=="__main__":
	main()
