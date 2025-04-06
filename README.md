# Emergence of novel SARS-CoV-2 variants in the Netherlands
This is the repository to accompany our publication, "Emergence of novel SARS-CoV-2 variants in the Netherlands" in *Scientific Reports*. In this repository, you can find the scripts we used to conduct a study of SARS-CoV-2 genomes to explore the viral population diversity in the Netherlands, within a global context. Our work lends insight into the genetic variation of SARS-CoV-2 in the later stages of the pandemic in April and early May.


## Dependencies
- mafft
- biopython
- numpy
- seaborn
- iqtree

## Usage
### Install dependencies using Conda
1. Install Anaconda or miniconda (if not already present on your system)
2. Create a new environment based on the `environment.yml` file
```bash
conda env create -f environment.yml
```
3. Activate the environment:
```bash
source activate sarscov2
```

### Data retrieval and preprocessing, and multiple sequence alignment
1. Complete, high quality (number of undetermined bases less than 1% of the whole sequence) genome sequences of SARS-COV-2 that were isolated from human hosts only were obtained from GISAID, NCBI and China’s National Genomics Data Center (NGDC) on June 13th. The dataset contained 29,503 sequences with unique identifiers in total, including the Wuhan-Hu-1 reference sequence (accession ID NC_045512.2). The “Collection date” field was also extracted for all sequences, and it is referred to as “date” throughout this work. The acknowledgement table for GISAID sequences can be found in Supplementary file 2 and the full list of sequence identifiers for NCBI and NGDC records are provided in Supplementary file 3.

You can find the final collection of SARS-CoV-2 genomes (after preprocessing and cleaning the metadata) on 4TU Research Database here. The `sarscov2_seq.fasta` file contains DNA sequences of 29,503 SARS-CoV-2 genomes in fasta format and the `sarscov2_metadata.tsv` contains the corresponding metadata. It has the fields ...

2. All sequences were aligned against the Wuhan-Hu-1 reference using MAFFT (v7.46) with the FFT-NS-fragment option, and the alignment was filtered to remove identical sequences to obtain 24,365 non-redundant genomes.
```bash
mafft --retree 1 --thread 4 sarscov2_seq.fasta > sarscov2_mafft_aln.fasta
```

### SARS-CoV-2 descriptive statistics
1. To generate the figures that describe the total SARS-CoV-2 population as well as the distribution of clades, run the python script `plot-descriptive.py` 
```bash
python plot-descriptive.py data/sarscov2_metadata.tsv -o img
```

2. Plot mutation frequencies
```bash
python plot-mutationfreq.py data/sarscov2_metadata.tsv data/sarscov2_mutations_coronapp.tsv -o img
```

3. Plot the locations of mutations on the S and N protein sequences, respectively usind the script `plot-sn-mutations.py`. In addition to the metadata file and the coronapp mutations, you need the S protein sequence in your `data` folder as `data/s-prot_nuc.fasta` and the mutation annotations in a tabular file as `data/mutations_annotated.tsv`
```bash
python plot-sn-mutations.py data/sarscov2_metadata.tsv data/sarscov2_mutations_coronapp.tsv -o img
```
