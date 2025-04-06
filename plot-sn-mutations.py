"""
This script generates the figure 4 from our publication
It needs to location of SARS-CoV-2 metadata file, mutation file (obtained from coronapp)
and the output directory for images
"""
import os
import argparse
import warnings
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description="Plot the locations of mutations on the S and N "
                                     "protein sequences")
    parser.add_argument('metadata', help="SARS-CoV-2 metadata file")
    parser.add_argument('mutation_file', help="SARS-CoV-2 mutations from coronapp")
    parser.add_argument('-o', '--output-dir', type=Path, required=True,
                        help="Output directory to save the images.")

    args = parser.parse_args()
    output_dir = args.output_dir
    metadata_file = args.metadata
    mutation_file = args.mutation_file
    sprot_fasta = "data/s-prot_nuc.fasta"
    mut_annot_file = "data/mutations_annotated.tsv"
    os.makedirs(output_dir, exist_ok=True)
    
    # Load metadata file describing each SARS-CoV-2 genome
    metadf = pd.read_csv(metadata_file, index_col=0, header=0, sep='\t')
    mutdf = pd.read_csv(mutation_file, sep='\t', index_col=0, header=0)
    mut_annotdf = pd.read_csv(mut_annot_file, sep='\t', index_col=0, header=0)
    sprotnuc_seq_dna = SeqIO.read(sprot_fasta, 'fasta')

    plt.style.use('seaborn-muted')
    warnings.simplefilter(action='ignore')

    # Figure 4  
    print("Plotting Figure 4") 
    metadf.loc[:,'acc'] = metadf.index # To keep the accesion IDs as a separate column
    # We plot only the samples that have the full date available
    idx_ok = metadf[metadf.date.apply(lambda x: len(x)==10)].index.intersection(set(mutdf.index)) 
    idx_nl = metadf.loc[idx_ok][metadf.loc[idx_ok].country=='Netherlands'].index
    mut_freqdf = mutdf.loc[idx_nl].groupby(['protein', 'varclass']).refpos.value_counts()

    fig,ax = plt.subplots(2, 1, figsize=(12, 6))
    pos = sorted(mut_freqdf.loc(axis=0)['S', :, :].index.get_level_values('refpos').unique())
    ax[0].bar(pos - mut_annotdf.loc['S','start'], mut_freqdf.loc(axis=0)['S', :, pos].sum(level='refpos').loc[pos],
              color='k', width=6, log=True)
    ax[0].text(23415 - mut_annotdf.loc['S', 'start'], mut_freqdf.loc(axis=0)['S', :, 23403].sum(level='refpos'),
               s='A23403G\nD614G', horizontalalignment='left', verticalalignment='top')
    ax[0].set_ylabel('Number of mutations')
    seq_len = len(sprotnuc_seq_dna)
    ax[0].set_xlim(0, seq_len)
    ax[0].set_xticks([*range(0, seq_len, 60)])
    ax[0].set_xticklabels([*range(0, seq_len, 20)])
    ax[0].xaxis.set_tick_params(labelrotation=90)
    ax[0].set_title('Spike protein')
    ax[0].grid('on', alpha=.5)
    
    pos = sorted(mut_freqdf.loc(axis=0)['N', :, :].index.get_level_values('refpos').unique())
    ax[1].bar(pos - mut_annotdf.loc['N', 'start'], mut_freqdf.loc(axis=0)['N',:,pos].sum(level='refpos').loc[pos], 
              color='k', width=2.5, log=True)
    ax[1].text(28887 - mut_annotdf.loc['N', 'start'], mut_freqdf.loc(axis=0)['N', :, 28881].sum(level='refpos'), 
               s='RG203KR', horizontalalignment='left', verticalalignment='top')
    ax[1].set_xlabel('Aminoacid position')
    ax[1].set_ylabel('Number of mutations')
    seq_len = mut_annotdf.loc['N', ['start','end']].diff()[1]
    ax[1].set_xlim(0, seq_len + 2)
    ax[1].set_xticks([*range(0, seq_len + 21, 60)])
    ax[1].set_xticklabels([*range(0, seq_len + 21, 20)])
    ax[1].xaxis.set_tick_params(labelrotation=90)
    ax[1].set_title('Nucleocapsid protein')
    ax[1].grid('on', alpha=.5)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig4.png'), dpi=600, format='png')
    plt.savefig(os.path.join(output_dir, 'fig4.pdf'), dpi=600, format='pdf')

if __name__ == "__main__":
    main()