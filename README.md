# Cell-type Specific Gene Networks and Drivers in Rheumatoid Arthritis

<img align="right" src="https://github.com/Aurelien-Pelissier/RA-drug-discovery/blob/main/img/PANDA.png" width=400>

The increasing number of available large RNA-seq datasets, combined with genome-wide association studies (GWAS), differential gene expression (DEG) studies, and gene regulatory networks (GRN) analyses have to potential to lead to the discovery of novel therapeutics. Yet, despite this progress, our ability to translate GWAS and DEG analyses into an improved mechanistic understanding of many diseases remains limited, as different analyses often disregard information about the cell-types mediating the disease.

In this repository, We compute sample-specific GRNs which enables the use of statistical techniques to compare network properties between phenotypic groups. Then, we leverage this collection of networks to rank transcription factors (TFs) according to their contribution to the observed differential gene expression between RA and control. 

### This repository support our publications:

[[1]]([https://www.biorxiv.org/content/10.1101/2023.12.28.573505](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1428773/full)) Pelissier A, Laragione T, Martinez MR, & Gulko PS. Cell-specific gene networks and drivers in rheumatoid arthritis synovial tissues. Frontier in immunology (2024): 2023-12.

[[2]]([https://www.biorxiv.org/content/10.1101/2023.12.28.573506](https://www.life-science-alliance.org/content/8/1/e202402808.abstract)) Pelissier A*, Laragione T*, Harris C, Gulko PS & , Martinez MR. BACH1 as a key driver in rheumatoid arthritis fibroblast-like synoviocytes identified through gene network analysis. Life Science Alliance (2025): 2023-12.

[//]: <> (Pelissier A*, Laragione T*, Martinez MR, & Gulko PS. BACH1 as key regulator in RA 2023. Planned.)

&nbsp;

## Constructing cell-type specific gene regulatory network
In this work, Gene regulatory networks are bipartite graphs, with edges connecting TF and their target gene (TG). Each edge has a weight representing the probability of a regulatory interaction between the connected nodes.
Briefly, we use PANDA [3] that integrates gene expression data with prior knowledge about TF-binding motif and protein-protein interactions by optimizing the weights of edges in the networks with iterative steps. PANDA's input consists of:
- A gene expression matrix. In this study we used bulk RNA expression, available in `Data/RA_gene_expression`
- Prior knowledge about TF-TF interactions and TF binding motif, available in `Data/PANDA_prior_knowledge`
Run the script `src/PANDA_network.py` to compute the network and analyse their edges.



Applied to our data, PANDA produced fully connected and directed networks of TFs to their target genes, comprising 644 TFs and 18992 genes. Then, we used LIONESS [4] to estimate an individual gene regulatory network for each sample in the population, which we utilized to make differential analysis of their edges and identify key TF regulators.

<p align="center">
  <img src="https://github.com/Aurelien-Pelissier/RA-drug-discovery/blob/main/img/LIONESS.png" width=500>
</p>


&nbsp;

## Key driver analysis
We use mergeomics [5]. You need:
Some network for your analysis. We used the GIANT network, downloaded at https://giant-v2.princeton.edu/download/.
- Run the script `src/KDA_analysis.py`

&nbsp;

## Experimental validation
In our article we focused on Synovial fibroblast and detected FOSL1, THBS1 and CFH as potential novel key regulators.
We performed silencing experiment on RA cell line and provide the results in `Data/silencing_experiments/silencing_data.xlsx`
- Run the script `src/silencing_analysis.py` to run the statistical test and combine the p-values with the Brown-Fisher method.

&nbsp;

## Reference
[3] Glass, Kimberly, et al. "Passing messages between biological networks to refine predicted interactions." PloS one 8.5 (2013): e64832.

[4] Kuijjer, Marieke Lydia, et al. "Estimating sample-specific regulatory networks." Iscience 14 (2019): 226-240.

[5] Shu, Le, et al. "Mergeomics: multidimensional data integration to identify pathogenic perturbations to biological systems." BMC genomics 17.1 (2016): 1-16.


