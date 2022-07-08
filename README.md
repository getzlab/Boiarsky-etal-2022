# Boiarsky-etal-2022
Code to reproduce methods &amp; results from Boiarsky et. al., *Nature Communications* 2022

Ig_genes.ipynb details how we determined which genes fall in the immunoglobulin loci, to remove them from downstream analyses.

4a_puritywork-published.ipynb contains our analysis of sample purity (% tumor cells in sample) using our Bayesian purity model, as described in our method. 

4d_limma.ipynb contains our limma-voom differential expression analysis comparing malignant or pre-malignant pseudobulk samples vs. normal pseudobulk samples.

helper_functions_published.py contains functions that are used throughout the other notebooks included in the repo.
