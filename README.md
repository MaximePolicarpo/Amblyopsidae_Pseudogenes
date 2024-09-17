# Amblyopsidae_Pseudogenes

## I - Analyse of loss-of-function mutation in vision genes

All the analysis can be reproduced using the R script "Amblyopsidae_Pseudogenes.R". All necessary files are included in this github repository

- AMAS_concatenated_alignment.prot.aln.treefile : Species tree obtained with BUSCO genes and ASTRAL-III (based on maximum likelihood trees from CODON alignments)
- AMAS_concatenated_alignment.cds.aln.treefile : Species tree obtained with BUSCO genes and ASTRAL-III (based on maximum likelihood trees from PROTEIN alignments)
- Stop_codon_table.tsv : Table of premature stop codon detected in Percopsiformes vision genes, with their position on the coding sequence
- Frameshift_table.tsv : Table of frameshits detected in Percopsiformes vision genes, with their nature (deletion/insertion), their size and their position on the coding sequence
- OtherLoF_table.tsv : Table of poroper start losses, proper stop losses, and splice site mutations detected in Percopsiformes vision genes
- IGV_verif_table.tsv : Table listing, for each found LoF mutation and for each individual, the number of reads supporting or not the mutation. 
- Species_table.tsv : Table of species investigated and their BUSCO results 
- Vision_genes_dNdS.csv  : Table of dN/dS for each vision genes and each species.
- Vision_genes_dNdS.Gadus.csv : Table of dN/dS for each vision genes and each species, which were computed with G. morhua sequences included in the alignments
- Vision_genes_dNdS.IND.csv : Table of dN/dS for each vision genes and each individual
- Vision_sequences_table.tsv : Table of vision gene sequences, their position and their state (Complete, truncated or pseudogene)
- Chase_CAVE_NODE_123SUM.tre : Species tree obtained with ultraconserved elements
  
