# CoExpPhylo
Collection of scripts associated with an integrated co-expression and phylogenetic analysis approach. The script can continue based on exisiting files. 


## Combination of co-expression and orthology

Based on defined genes in different species, this script searches for co-expressed genes in each species. The predefined genes of different species should be orthologs. To find out if the co-expressed genes are also orthologs, a phylogenetic tree is constructed. This allows to harness the available transcriptome data sets across species borders. If gene expression networks are conserved across species, the sequences identified through co-expression should cluster in phylogenetic trees that are constructed in the final step.

**Please note that the script can continue based on exisiting files. If you want to start a new analysis, the output-folder should be empty or not exist yet.**

## Installation via Conda environment

1. Clone the repository:
```
git clone https://github.com/bpucker/CoExpPhylo.git
cd CoExpPhylo
```
2. Create and activate the Conda environment
```
conda env create -f environment.yml
conda activate env_CoExpPhylo
```

### Additional Requirement for automatic upload to iTOL ([Letunic and Bork (2024)](https://doi.org/10.1093/nar/gkae268))
If you want to use the upload to iTOL feature, download the `iTOL_uploader.pl` script from: https://itol.embl.de/help.cgi#batch

Save it in your project directory and ensure it is executable:
```
chmod +x iTOL_uploader.pl
```

Eventually, additional perl packages must be installed to execute the `iTOL_uploader.pl` script.


## Usage

```
python3 coexp_phylo.py --config <FILE> --out <DIR> --mode <MODE>

Mandatory arguments:
  --config             STR    Config file.
  --out                STR    Directory for temporary and output files.
 
		
Optional arguments:
 ----
  ANNOTATION
 ----
  --anno               STR     Annotation file matching reference
  --reference          STR     Annotation peptide file
  --seqs_cluster_anno  FLOAT   Percentage of sequences per cluster used for annotation [50.0]

 ----
  PER-SPECIES COEXPRESSION ANALYSIS
 ----
  --r                  FLOAT   Correlation coefficient cutoff [0.7]
  --p                  FLOAT   P-value cutoff [0.05]
  --numcut             INT     Number of co-expressed genes [100]
  --min_exp_cutoff     INT     Minimal expression per gene[30]
  --mindetect          INT     Minimal number of coexpressed baits for a seqeunce to be included in the anaylsis[1]

 ----
  DIAMOND
 ----
  --cpub               INT     Number of cores to use [4]
  --batch              INT     Number of Blastp jobs to run in parallel [7]
  --scorecut           FLOAT   Minimal BLAST hit score cutoff [100.0]
  --simcut             FLOAT   Minimal BLAST hit similarity cutoff [80.0]
  --lencut             INT     Minimal BLAST hit length cutoff [100]
  --evalue             FLOAT   E-value cutoff [1e-5]
  
 ----
  CLUSTERING
 ----
  --minseqcutoff       INT     Minimal sequence cluster size [10]
  --mincoexpseqcutoff  INT     Minimal number of species with coexpressed genes to form cluster [3]

 ----
  ALIGNMENT
 ----
  --alnmethod          STR     Alignment algorithm (mafft|muscle)[mafft] 
  --mafft              STR     Full path to MAFFT [mafft]
  --muscle             STR     Full path to muscle [muscle]
  --occupancy          FLOAT   Minimal occupancy for a position during alignment to be kept [0.1]

 ----
  TREE CONSTRUCTION
 ----
  --treemethod         STR     Tree construction algorithm (fasttree|raxml|iqtree) [fasttree]
  --raxml              STR     Full path to RAxML [raxml]
  --fasttree           STR     Full path to FastTree [fasttree]
  --iqtree             STR     Full path to IQ-TREE [iqtree]
  --cpur               INT     Number of cores for tree construction [cpub]
  --clean              BOL     If argument is given: Cleaning of ids for correct visalization in iTOL

 ----
  BATCH UPLOAD ITOL
 ----
 if batch upload is wanted, the following arguments are mandatory
  --API <API-KEY>
  --proj_name <PROJECT_NAME>(project must already exist in your iTOL account)
  optional:
  --upload_script <PATH_TO_iTOL_UPLOADER.PL>

```

### Mandatory arguments
`--config` specifies a config file that contains all the information about the input files. The columns in this file need to be comma-separated. Each row describes one data set. The columns are: ID, TPM file, CDS file, name of file with IDs of bait sequences, and optional PEP file.

ID = This is usually the species name, but should not contain any spaces or other weired characters. Using a-Z and 0-9 is fine with underscores as space replacements. This ID is used as a prefix for all sequences of this species to avoid ambiguities when combining the sequences of different species in one file.

TPM file = First column contains the gene IDs and the first row contains the sample names. All other fields in this table are gene expression values. The IDs in this file need to match the sequences in the CDS file and also the IDs given as baits.

CDS file = This is a multiple FASTA file with the coding sequences of this species. The sequence IDs need to match the TPM file and also the IDs specified as baits in the fourth input file type.

Bait sequence ID file = This file contains the ID of genes that should be used as baits. For example, these could be known genes involved in the upstream part of a pathway if the objective is to discover genes further downstream in the pathway. One ID should be given per line.

PEP file = This is a multiple FASTA file with the peptide sequences of this species. The sequence IDs need to match the TPM file and also the IDs specified as baits in the fourth input file type. The specification of a PEP file is optional, but speed up the analyses. Hence, the config file can also contain four columns.


`--out` specifies the output folder. All temporary file and the final output files will be placed in this folder. The folder will be created if it does not exist already.


### Optional arguments

#### Annotation
`--anno` specifies a tab-separated annotation file. IDs need to be located in the first column and the annotation text need to be located in the second column.

`--reference` specifies the peptide sequence file as reference for the analysis. This allows an effective annotation in the final steps.

`--seqs_cluster_anno` specifies the percentage of sequences per cluster that should be used for annotation. This allows an improvement in time.

#### Per-species coexpression analysis
`--r` specifies the minimal correlation coefficient for genes to be considered. The default value is 0.7.

`--p` specifies the maximal p-value in the correlation calculation for genes to be considered. The default value is 0.05.

`--numcut` specifies the maximal number of genes to consider in the co-expression analysis. Only these top sequences are analyzed in the next steps. The default value is 100. An increase of this number will substantially increase the run time.

`--min_exp_cutoff` specifies the minimal combined expression of a any given gene to be considered in the co-expression analysis. Default: 30.

`--mindetect` specifies the minimal number of baits that any given gene needs to be co-expressed with. Depending on the number of baits, this value can be increased. A higher value will generally boost the specificity at the cost of sensitivity. The strictes possible option is setting this number to the number of baits. Default: 1.

#### DIAMOND
`--cpub` specifies the number of cores for DIAMOND blastp and other operations. The default value is 4.

`--batch` specifies the number of blastp jobs to run in parallel. The default value is 7.

`--scorecut` specifies the minimal BLAST score cutoff for hits to be considered. The default value is 100. The value is relevant for the construction of sequence clusters once co-expressed genes have been identified.

`--simcut` specifies the minimal similarity of a DIAMOND blast hit to be considered. The default value is 80.

`--lencut` specifies the minimal lenght of a DIAMOND blast hit to be considered. The default value is 100.

`--evalue` specifies the maximal e-value of a DIAMOND blast hit to be considered. The default value is 10<sup>-5</sup>.

#### Clustering
`--minseqcutoff` specifies the minimal number of sequences that must belong to a cluster to trigger the construction of a phylogenetic tree. Extremely small clusters are masked to avoids inflation of cluster numbers. Default: 10.

`--mincoexpseqcutoff` specifies the minimal number of species with co-expressed genes that must be in the same sequence cluster (tree). Default: 3.

#### Alignment
`--alnmethod` specifies the algorithm/tool used for construction of the multiple sequence alignment as basis of the phylogenetic tree construction. Currently, MAFFT and MUSCLE are the supported options. Default is MAFFT.

`--mafft` specifies the MAFFT path. This option can be used if MAFFT is not in the `$PATH` variable or if a specific version should be used. The default is 'mafft'.

`--muscle` specifies the MUSCLE path. This option can be used if MUSCLE is not in the `$PATH` variable or if a specific version should be used. The default is 'muscle'.

`--occupancy` specifies the minimum column occupancy threshold for the multiple sequence alignment. Only alignment columns where at least the specified fraction of sequences contain a residue (i.e., are not gaps `-`) are retained for further analysis. The default value is 0.1.

#### Tree construction
`--treemethod` specifies the algorithm/tool used for construction of the phylogenetic trees once sequence clusters are identified. Currently, FastTree, RAxML, and IQ-TREE are the supported options. Default is FastTree.

`--raxml` specifies the RAxML path. This option can be used if RAxML is not in the `$PATH` variable or if a specific version should be used. The default is 'raxml'.

`--fasttree` specifies the full path to FastTree. This option can be used if FastTree is not in the `$PATH` variable or if a specific version should be used. The default is 'fasttree'.

`--iqtree` specifies the full path to IQ-TREE. This option can be used if IQ-TREE is not in the `$PATH` variable or if a specific version should be used. The default is 'iqtree'.


#### Batch uplaod to iTOL ([Letunic and Bork (2024)](https://doi.org/10.1093/nar/gkae268))
If wanted, the trees can automatically be uploaded to iTOL. To use this option, you must have an active standard subscription. 

##### Mandatory arguments
`--API` specifies the API key.
`--pro_name` specifies the project name. The project name should be unique among your workspaces. The project must already exist.
##### Optional arguments
`--upload_script` specifies the path to the uplaod perl script. Per default, the script is downloaded included in this repository and a provision of the path is not necessary. The script can be downloaded via [itol.embl.de](https://itol.embl.de/help.cgi). If you downloaded the script manually, you might need to provide the path.


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

## Computational requirements

The pipeline performs several computationally intensive steps, including a coexpression analysis, all-vs-all sequence comparisons, and multiple sequence alignments.

Concise requirements depend on the selected tools (e.g., MAFFT vs. MUSCLE), the number of input species, and size and complexity of the input datasets. Below, we provide an overview of the hardware used and corresponding runtimes for different datasets:

|Data|CPU cores|RAM|Runtime (default parameters)|
|---|---:|---:|---|
|Example dataset (small test case)|28|256 GB|~1 minute|
|Dataset as described in publication|28|256 GB|~5 hours|

These examples were run on a high-performance computing cluster. For smaller datasets or fewer species, lower resources (e.g., 16-32 GB RAM and 8-16 cores) may suffice, though runtimes will increase accordingly.

## License

GNU GENERAL PUBLIC LICENSE, Version 3

## References

This repository.

Grünig, N. & Pucker, B. (2025). CoExpPhylo – A Novel Pipeline for Biosynthesis Gene Discovery. bioRxiv; doi: [10.1101/2025.04.03.647051](https://doi.org/10.1101/2025.04.03.647051).

Ivica Letunic, Peer Bork, Interactive Tree of Life (iTOL) v6: recent updates to the phylogenetic tree display and annotation tool, Nucleic Acids Research, Volume 52, Issue W1, 5 July 2024, Pages W78–W82, https://doi.org/10.1093/nar/gkae268
