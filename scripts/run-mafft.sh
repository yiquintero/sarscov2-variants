#!/bin/sh

#SBATCH --qos=short
#SBATCH --partition=general
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4GB
#SBATCH --mail-type=END

srun mafft --retree 1 --thread 4 data/sarscov2_sequences.fasta > output/mafft-aln.fasta

