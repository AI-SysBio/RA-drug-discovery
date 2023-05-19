The lists in this folder are potential key drivers TF and target genes, among with their expression level in each cell type (from 1rst to 10th quantile). It is the result from the combination of two independent pipelines. Depending on public availability, not all cell types could be studies by both methods. Therefore, the list of key driver genes obtained by each cell type are defined differently depending on the method used.

    - Method 1 relies on previously published cell type gene regulatory networks and look for specific genes associated to many RA-associated signatures. The set of cell type specific key drivers is defined as genes that are found in the KDA list from both the DEG list and Litterature list.
    
    - Method 2 build cell type specific gene regulatory network in synovial fibroblast with B cells, T cells, Fibroblasts and monocytes. Since these networks are bipartite graph between TFs and their respective target genes, the list is build differently for TF and their target genes. Key driver TFs are defined by looking at the TFs maximizing the correlation between the differential expression of their edge weights to their target with the gene differential expression in the cell type. While not much information can be extracted from key driver target genes in this context, we contrain the set of key driver target genes (to minimize false positive) with the condition that they should have a t-statistic between RA and control gene expression > 1.
	
	- Method 1 + Method 2 combine the two pipelines above and keeps the overlap of the two as key drivers.
	
	- For cell types where method 1 is not available, we can still contrain the list to genes found as key drivers in many other cell types.
	  We chose genes that are found as key driver consistently in other cell types to reduce false positive.

  Below is a table of available method for each cell type. Obviously, Lists of cell types where both methods can be combined can be considered more reliable with less false positive

  Celltype 		Method1		Method2
  ------------------------------------
  Adipocyte			YES			NO
  B cell			YES			YES
  DC				YES			NO
  Eosinophil		YES			NO
  Fibroblast		NO			YES
  Macrophage		YES			NO
  Monocyte			YES			YES
  NK cell			YES			NO
  T cell			YES			YES

  See the list of networks for method 1 here (https://hb.flatironinstitute.org/download)