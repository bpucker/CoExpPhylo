# CoExpPhylo
Collection of scripts associated with co-expression and phylogenetic analysis approach



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


