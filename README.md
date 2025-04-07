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

### Data retrieval and preprocessing
1. Complete, high quality (number of undetermined bases less than 1% of the whole sequence) genome sequences of SARS-COV-2 that were isolated from human hosts only were obtained from GISAID, NCBI and China’s National Genomics Data Center (NGDC) on June 13th. The dataset contained 29,503 sequences with unique identifiers in total, including the Wuhan-Hu-1 reference sequence (accession ID NC_045512.2). The “Collection date” field was also extracted for all sequences, and it is referred to as “date” throughout this work. The acknowledgement table for GISAID sequences can be found in Supplementary file 2 and the full list of sequence identifiers for NCBI and NGDC records are provided in Supplementary file 3.

You can find the final collection of SARS-CoV-2 genome sequences (after preprocessing and cleaning the metadata) on 4TU Research Database [here](https://doi.org/10.4121/11bff1ea-4784-463e-90d0-eb2e2b64fe96). The `sarscov2_sequences.fasta` file contains DNA sequences of 29,503 SARS-CoV-2 genomes in fasta format. In addition, we provide the genome metadata file (`sarscov2_metadata.tsv`), mutations obtained from the coronapp web application (`sarscov2_mutations_coronapp.tsv`) on the same 4TU.ResearchData database. On this GitHub repository, you can also find the annotations for the mutations (`data/mutations_annotated.tsv`) and the nucleotide sequence for the S protein in fasta format (`data/s-prot_nuc.fasta`).

At the start, you should have all the files below under the same `data` folder:

- `sarscov2_metadata.tsv`: This is tab-separated-file that contains all the metadata related to the 29,503 SARS-CoV-2 genomes used in our work. The table is indexed by the genome sequence IDs and consists of 27 fields. 
- `sarscov2_sequences.fasta`: DNA sequences of 29,503 SARS-CoV-2 genomes in fasta format.
- `sarscov2_mutations_coranapp.tsv`: After running the `mafft` alignment on all SARS-CoV-2 sequences, we uploaded the alignment file after we filtered to remove the identical sequences and trimmed to remove the gaps from the Wuhan-Hu-1 reference (accession ID NC_045512.2) to the *coronapp* web application. Since the *coronapp* application has recently enforced a limit on maximum number of sequences allowed, we provide the original output in the 4TU.ResearchDatabase.
- `mutations_annotated.tsv`: This is a tab-separated file that lists all the mutations that had annotations avaialable at the time of our analyses. It lists each mutation's location, locustag, ID of the protein it is located in and the result of the mutation.
- `s-prot_nuc.fasta`: The nucleotide sequence of the S protein in fasta format.

### Plotting the results from variant analysis
1. To generate Figures 1 and 2 from our publication use the python script `plot-descriptive.py`: these figures describe the total SARS-CoV-2 population per continent as well as the distribution of clades over time in the Netherlands, the UK, Australia, Singapore and China.
```bash
python plot-descriptive.py data/sarscov2_metadata.tsv -o img
```

2. To generate Figure 3 from our publication use the python script `plot-mutationfreq.py`: this figure shows how the number of mutations per sample per day change over time in the Netherlands, Europe and globally.
```bash
python plot-mutationfreq.py data/sarscov2_metadata.tsv data/sarscov2_mutations_coronapp.tsv -o img
```

3. To generate Figure 4 use the python script `plot-sn-mutations.py`: this figure shows the total number of mutations in the S and N proteins in samples from the Netherlands. In addition to the metadata file and the coronapp mutations, you need the S protein sequence in your `data` folder as `data/s-prot_nuc.fasta` and the mutation annotations in a tabular file as `data/mutations_annotated.tsv`
```bash
python plot-sn-mutations.py data/sarscov2_metadata.tsv data/sarscov2_mutations_coronapp.tsv -o img
```

4. To generate Figures 5 and 6 use the python script `plot-global-mutations.py`: Figure 5 shows the total frequency of the top 15 mutations observed in the most sampled countries in our dataset in total, while Figure 6 shows how the frequency of the same mutations change over time.
```bash
python plot-global-mutations.py data/sarscov2_metadata.tsv data/sarscov2_mutations_coronapp.tsv -o img
```

5. To generate Figure 7 use the python script `plot-mutationfreq-time.py`: Figure 7 shows the change in frequency of the top 15 mutations in the Netherlands averaged over a period of 7 days. The top 15 mutations are the same as those in Figures 5 and 6 to be consistent.
```bash
python plot-mutationfreq-time.py data/sarscov2_metadata.tsv data/sarscov2_mutations_coronapp.tsv -o img
```

### Multiple sequence alignment and phylogenetic tree analysis

1. Align all sequences against the Wuhan-Hu-1 reference using MAFFT (v7.46) with the FFT-NS-fragment option
```bash
mafft --retree 1 --thread 4 data/sarscov2_seq.fasta > output/sarscov2_mafft_aln.fasta
```
- We recommend submitting the `mafft` alignment as a job to a HPC, you can find an example script for SLURM jobs in `scripts/run-mafft.sh`
2. Filter the alignment to remove gaps and identical sequences to obtain 24,365 non-redundant genomes using `scripts/alignment_utils.py`. This will produce an output file titled `output/sarscov2_mafft_aln_refined_trimmed.fasta`
```bash
python scripts/alignment_utils.py output/sarscov2_mafft_aln.fasta
```
3. Generate a phylogenetic tree using IQTREE To generate a tree, use IQTREE (v2.05) with GTR model, allowing to collapse non-zero branches, and ultrafast bootstrap with 1000 replicates:
```bash
iqtree -czb -m GTR -nt AUTO -ntmax 4 -pre output/iqtree-tree_collapsed -s output/sarscov2_mafft_aln_refined_trimmed.fasta -seed 1 -v
```
4. Finally, you can use [iTOL](https://itol.embl.de/) to draw the phylogenetic tree.


