# inferCNVpy_pipeline

## Intro

This project uses inferCNVpy tool to conduct CNV analysis per sample basis from a single AnnData object. 
The code snippets here are for doing the following:

* Preprocess AnnData object, annotate chromosome locations etc for the tool to run properly
* Run basic inferCNVpy analysis for each sample in your object and add the results in the original object
* Make standard plots like heatmaps and UMAPs
* Run linear association analysis for gene expression and cnv_score
* 


## To Do

- [ ] construct main script
- [ ] add plots from each part of the analysis
- [ ] write walkthrough for the pipeline using sample data
- [ ] Make the optional use aggregated controlset and opt the cells out
- [ ] fork the original
