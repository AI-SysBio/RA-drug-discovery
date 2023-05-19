These were downloaded from the GRAND database
https://grand.networkmedicine.org/

Example: (see PPI and regulation priors in the link below)
https://grand.networkmedicine.org/tissues/Spleen_tissue/
Note: the same prior data were used for all the tissues currently in the database

[ref] Sonawane, Abhijeet Rajendra, et al. "Understanding tissue-specific 
      gene regulation." Cell reports 21.4 (2017): 1077-1088.


Briefly, these data were build as follow:

Transcription Factor Motif Information:
	We downloaded Homo sapiens transcription factor motifs with direct/inferred
	evidence from the Catalog of Inferred Sequence Binding Preferences (http://
	cisbp.ccbr.utoronto.ca/, accessed July 7, 2015). For each transcription factor,
	we selected the motif with the highest information content, and we mapped its
	position weight matrix to the human genome (hg19) using the Find Individual
	Motif Occurrences program (Grant et al., 2011). We retained significant
	hits (p < 1e-5) that occurred within the promoter ([-750, +250] around the
	transcriptional start site) of Ensembl genes (GRCh37.p13; annotations from
	https://genome.ucsc.edu/cgi-bin/hgTables, accessed September 3, 2015).
	We intersected this map with the expression data, resulting in a set of canonical
	regulatory interactions from 644 transcription factors to 30,243 genes,
	which we used to construct our regulatory network models.

Prior Protein-Protein Interaction Data:
	A protein-protein interaction network between transcription factors in our motif
	prior was constructed based on interaction scores from StringDb version (v.)10
	(https://string-db.org, accessed October 27, 2015).