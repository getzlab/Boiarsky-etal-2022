# Boiarsky-etal-2022
Code to reproduce methods &amp; results from Boiarsky et. al., *Nature Communications* 2022

For now, please cite our preprint, available on medRxiv!
```
Single Cell Characterization of Myeloma and its Precursor Conditions Reveals Transcriptional Signatures of Early Tumorigenesis
Rebecca Boiarsky, Nicholas J. Haradhvala, Jean-Baptiste Alberge, Romanos Sklavenitis-Pistofidis, Tarek H Mouhieddine, 
Oksana Zavidij, Ming-Chieh Shih, Danielle Firer, Mendy Miller, Habib El-Khoury, Shankara K. Anand, François Aguet, 
David Sontag, Irene M. Ghobrial, Gad Getz
medRxiv 2022.02.01.22270128; doi: https://doi.org/10.1101/2022.02.01.22270128
```

## System Requirements
Our code was run using the following software:
- Python version 3.8.13
- R version 4.1.0

#### Python packages to install:
- scanpy==1.7.1
- datatable
- scikit_posthocs
- statsmodels
- scipy
- seaborn
- re

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
To create an anndata object from the raw data publicly available on GEO (accession number GSE193531) and the supp/source tables published with the paper, follow the code in the notebook [0_reproduce_results_from_raw_data.ipynb](https://github.com/getzlab/Boiarsky-etal-2022/blob/master/0_reproduce_results_from_raw_data.ipynb) (this notebook is being updated, further sections will be completed soon).

## Analysis notebooks:
The notebook 0_reproduce_results_from_raw_data.ipynb is a good starting place if one wants to explore our data themselves. More detailed versions of our analysis code is contained in the following notebooks:
- Ig_genes.ipynb details how we determined which genes fall in the immunoglobulin loci, to remove them from downstream analyses.
- 4a_puritywork-published.ipynb contains our analysis of sample purity (% tumor cells in sample) using our Bayesian purity model. It also contains the code to generate Fig. 2a from our paper. 
- 4d_limma.ipynb contains our limma-voom differential expression analysis comparing malignant or pre-malignant pseudobulk samples vs. normal pseudobulk samples.
- 5_NMF_rawdata-moreHVG-published.ipynb contains code related to generating the input data for SignatureAnalyzer and our analysis of SignatureAnalyzer results, including Figs. 3a-d from our paper.
- 5b_heterogeneity.ipynb contains code for analyzing the heterogeneity of signature expression within tumor samples, including Fig. 4c from our paper.
- helper_functions_published.py contains functions that are used throughout the other notebooks included in the repo.
