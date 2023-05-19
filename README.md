# Cell-type Specific Gene Networks and Drivers in Rheumatoid Arthritis

<img align="right" src="https://github.com/Aurelien-Pelissier/RA-drug-discovery/blob/main/img/PANDA.png" width=400>

The increasing number of available large RNA-seq datasets, combined with genome-wide association studies (GWAS), differential gene expression (DEG) studies, and gene regulatory networks (GRN) analyses have to potential to lead to the discovery of novel therapeutics. Yet, despite this progress, our ability to translate GWAS and DEG analyses into an improved mechanistic understanding of many diseases remains limited, as different analyses often disregard information about the cell-types mediating the disease.

In this repository, We compute sample-specific GRNs which enables the use of statistical techniques to compare network properties between phenotypic groups. Then, we leverage this collection of networks to rank transcription factors (TFs) according to their contribution to the observed differential gene expression between RA and control. 

### This repository support our publication:

[1] Pelissier A*, Laragione T*, Martinez MR, & Gulko PS. Cell-type Specific Gene Networks and Drivers in Rheumatoid Arthritis (2023). Planned.

[//]: <> (Pelissier A*, Laragione T*, Martinez MR, & Gulko PS. BACH1 as key regulator in RA 2023. Planned.)


## Constructing cell-type specific gene regulatory network
We use PANDA and LIONESS:
- Cell-specific gene expression matrix. In this study we used bulk but you can also use single cell and average them out to make it equivalent, available in `Data/RA_gene_expression`
- Prior knowledge about TF-TF interactions and TF binding motif, available in `Data/PANDA_prior_knowledge`
Run the script `src/PANDA_network.py` to compute the network and analyse their edges.

<p align="center">
  <img src="https://github.com/Aurelien-Pelissier/RA-drug-discovery/blob/main/img/LIONESS.png" width=300>
</p>

## Key driver analysis
We use mergeomics src/KDA. You need:
Some network for your analysis. We used the GIANT network, downloaded at https://giant-v2.princeton.edu/download/.
- Run the script `src/KDA_analysis.py`

## Experimental validation
In our article we focused on Synovial fibroblast and detected FOSL1, THBS1 and CFH as potential novel key regulators.
We performed silencing experiment on RA cell line and provide the results in `Data/silencing_experiments/silencing_data.xlsx`
- Run the script `src/silencing_analysis.py` to run the statistical test and combine the p-values with the Brown-Fisher method.

## Reference
[2] Shu, Le, et al. "Mergeomics: multidimensional data integration to identify pathogenic perturbations to biological systems." BMC genomics 17.1 (2016): 1-16.

[3] Kuijjer, Marieke Lydia, et al. "Estimating sample-specific regulatory networks." Iscience 14 (2019): 226-240.
