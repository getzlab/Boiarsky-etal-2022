# Boiarsky-etal-2022
Code to reproduce methods &amp; results from Boiarsky et. al., *Nature Communications* 2022

## System Requirements
Our code was run using the following software:
- Python version 3.8.13
- R version 4.1.0

#### Python packages to install:

pip installs:
- scanpy==1.7.1

other:
- clone the [SignatureAnalyzer-GPU](https://github.com/broadinstitute/SignatureAnalyzer-GPU) repo to run Signature Analyzer NMF decomposition

#### R packages to install:
- stringr_1.4.0 
- tibble_3.1.7  
- Matrix_1.4-1  
- edgeR_3.36.0  
- ggplot2_3.3.6
- tidyr_1.2.0   
- dplyr_1.0.9   
- limma_3.50.3

## To create an anndata object containing our single cell RNAseq data:
Raw UMI counts & metadata are available on GEO, accession number GSE193531

## Notebooks:
- Ig_genes.ipynb details how we determined which genes fall in the immunoglobulin loci, to remove them from downstream analyses.
- 4a_puritywork-published.ipynb contains our analysis of sample purity (% tumor cells in sample) using our Bayesian purity model, as described in our method. 
- 4d_limma.ipynb contains our limma-voom differential expression analysis comparing malignant or pre-malignant pseudobulk samples vs. normal pseudobulk samples.
- helper_functions_published.py contains functions that are used throughout the other notebooks included in the repo.
