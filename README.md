# CoExpPhylo
Collection of scripts associated with co-expression and phylogenetic analysis approach.


## Combination of co-expression and orthology

Based on defined genes in different species, this script searches for co-expressed genes in that species. The predefined genes of different species should be orthologs. To find out if the co-expressed genes are also orthologs, a phylogenetic tree is constructed. This allows to harness the available transcriptome data sets across species borders. If gene expression networks are conserved across species, the sequences identified through co-expression should cluster in phylogenetic trees that are constructed in the final step.


```
Usage:
  python3 coexp_phylo.py --config <FILE> --out <DIR>

Mandatory:
  --config   STR    Config file.
  --out      STR    Directory for temporary and output files.
 
		
Optional:
  --anno     STR     Annotation file matching reference
  --araport  STR     Araport11 peptide file
  --r        FLOAT   Correlation coefficient cutoff [0.7]
  --p        FLOAT   P-value cutoff [0.05]
  --numcut   INT     Number of co-expressed genes [100]
  --cpu      INT     Number of cores to use [4]
  --scorecut FLOAT   Minimal BLAST hit score cutoff [100.0]
  --simcut   FLOAT   Minimal BLAST hit similarity cutoff [60.0]
  --lencut   INT     Minimal BLAST hit length cutoff [100]
  --mode     STR     Tree construction algorithm (fasttree|raxml) [fasttree]
  --mafft    STR     Full path to MAFFT [mafft]
  --raxml    STR     Full path to RAxML [raxml]
  --fasttree STR     Full path to FastTree [fasttree]
  --cpur     INT     Number of cores for tree construction [cpu]
```


`--config` specifies a config file that contains all the information about the input files. The columns in this file need to be comma-separated. Each row describes one data set. The columns are: ID, TPM file, CDS file,file with IDs of bait sequences.

ID = This is usually the species name, but should not contain any spaces or other weired characters. Using a-Z and 0-9 is fine with underscores are space replacements. This ID is used as an prefix for all sequences of this species to avoid ambiguities when combining the sequences of different species in one file.

TPM file = First column contains the gene IDs and the first row contains the sample names. All other fields in this table are gene expression values. The IDs in this file need to match the sequences in the CDS file and also the IDs given as baits.

CDS file = This is a multiple FASTA file with the coding sequences of this species. The sequence IDs need to match the TPM file and also the IDs specified as baits in the last input file type.

Bait sequence ID file = This file contains the ID of genes that should be used as baits. For example, these could be known genes involved in the upstream part of a pathway if the objective is to discover genes further downstream in the pathway. One ID should be given per line.








## Script for the annotation of sequence clusters

```
Usage:
  python3 annotate_seqs.py --in <FILE> --out <DIR> --ref <FILE>

Mandatory:
  --in   STR         Query multiple FASTA file. 
  --ref  STR         Reference multiple FASTA file
  --out  STR         Directory for temporary and output files.
 
		
Optional:
  --anno STR         Annotation file matching reference
```


`--in` specifies a multiple FASTA file with sequences that should have a similar function. An annotation will be assigned to this group of sequences.

`--ref` specifies a multiple FASTA file that contains the reference sequences. A functional annotation of these sequences is available and will be used to annotate the query sequences.

`--out` specifies the output folder. This folder will be created if it does not exist already.

`--anno` specifies the annotation file. This file contains a tab-delimited table with the reference sequence ID in the first column and the annotation in the following column(s). All following columns will be merged with a ";" as separator.


