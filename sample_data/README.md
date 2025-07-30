# Sample Dataset for CoExpPhylo

This folder contains a synthetic dataset created to demonstrate and test the functionality of **CoExpPhylo**. 

It does **not represent real biological data** and is intended for testing, debugging, and illustrating the required input data formats.
The dataset is structured to mimic a typical multi-species coexpression and phylogenetic analysis setup.

## Folder Structure and Contents
```css
example_dataset/
├── sample_config.csv
├── data/
│   ├── Species_A/
│   │   ├── Species_A.fasta
│   │   ├── Species_A_counts.tsv
│   │   └── bait.txt
│   ├── Species_B/
│   │   ├── Species_B.fasta
│   │   ├── Species_B_counts.tsv
│   │   └── bait.txt
│   └── ... (Species_C to Species_L)
└── reference/
    ├── reference.fasta
    ├── reference.pep
    └── reference_annotation.tsv
 
```
#### `sample_config`
A configuration file listing all species and paths of input data for CoExpPhylo. **This file is required to run CoExpPhylo.**

Before running the program, update all file paths in this configuration file to match your environment.

#### `data/` 
Contains one subfolder for each species (`Species_A` to `Sepcies_L`), with the following files:

- `<Species>.fasta` - Sytnethic transcript sequences of the species.
- `<Species>_counts.tsv`- Simulated gene expression count table:
    - Rows: Gene IDs
    - Columns: Sample IDs
- `bait.txt` - Species-specific ist of bait gene(s) used as starting points for coexpression analysis. 

#### `reference/`
Contains synthetic reference data to annotate output data:
- `reference.fasta`- Reference transcript sequences.
- `reference.pep` - Reference protein sequences (translated `reference.fasta`).
- `reference_annotation.tsv`- Tab-separated annotation table mapping reference genes to synthetic functional descriptions.


## Usage with CoExpPhylo
Please ensure that you **adjust the paths** in `sample_config.csv` before running the workflow.

```bash
python3 /path/to/coexp_phylo.py \
--config /path/to/sample_data/sample_config.csv \
--out path/to/output/  \
--mode batch \
--reference /path/to/sample_data/data/reference/reference.pep \
--anno /path/to/sample_data/data/reference/reference_annotation.tsv
```

## Output Files
This folder contains a collection of output files that would be generated using this dataset and CoExpPhylo with default parameters:
- `functional_annotation_of_clusters.txt` - Overview about the generated OCGs and their number of sequences as well as the annotation.
- `aln` - Folder containg all seqeunces belonging to the individual OCGs.
- `trees` - Folder containing all `.tree` files in Newick format. These files could be uploaded to iTOL.


## Data Origin
All files were synthetically generated and do not represent real organisms or biological measurements. They are meant solely for demonstration purposes and are free to use for testing.