#!/usr/bin/env python3
"""
This script generates the figures 3 and 4 from our publication
It needs to location of SARS-CoV-2 metadata file, mutation file (obtained from coronapp)
and the output directory for images
"""
import os
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from datetime import datetime

def main():
    parser = argparse.ArgumentParser(description="Plot mutation frequencies of SARS-CoV-2 genomes in the Netherlands over time"
                                     "and the locations of these mutaitons on the genomedescriptive statistics of the total "
                                     "SARS-CoV-2 population and the individual clades")
    parser.add_argument("metadata", help="SARS-CoV-2 metadata file")
    parser.add_argument('mutation_file', help="SARS-CoV-2 mutations from coronapp")
    parser.add_argument('-o', '--output-dir', type=Path, required=True,
                        help="Output directory to save the images.")

    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load metadata file describing each SARS-CoV-2 genome
    metadf = pd.read_csv(metadata, index_col=0, header=0, sep='\t')
    mutdf = pd.read_csv(mutation_file, sep='\t', index_col=0, header=0)

    plt.style.use('seaborn-muted')

    # Figure 3    
    cmap = {'green': '#66c2a5', 'orange': '#fc8d62', 'blue': '#8da0cb'}
    metadf.loc[:,'acc'] = metadf.index # To keep the accesion IDs as a separate column
    # We plot only the samples that have the full date available
    idx_ok = metadf[metadf.date.apply(lambda x: len(x)==10)].index.intersection(set(mutdf.index)) 
    idx_nl = metadf.loc[idx_ok][metadf.loc[idx_ok].country=='Netherlands'].index
    mut_freqdf = mutdf.loc[idx_ok].groupby(['sample','varclass']).count().refpos
    dates_axis = sorted(metadf.loc[idx_ok].date.unique())
    
    datedf = pd.concat([metadf.loc[idx_nl].date, mut_freqdf.loc(axis=0)[idx_nl, :].sum(level='sample')],axis=1)
    datedf = datedf.set_index('date').sum(level='date')
    add_dates = set(dates_axis).difference(datedf.index)
    add_dates = pd.DataFrame([np.nan]*len(add_dates), index=add_dates, columns=['refpos'])
    datedf = datedf.append(add_dates)
    datedf.sort_index(inplace=True)
    tot_sum = metadf.loc[idx_nl].groupby('date').count().acc.append(add_dates.refpos)
    solid = [date for (date, val) in zip(dates_axis, tot_sum.loc[dates_axis].values) if val>=5]
    faded = list(set(dates_axis).difference(solid))

    fig,ax = plt.subplots(figsize=(12, 4))
    ax.scatter([dates_axis.index(x) for x in solid], datedf.loc[solid].refpos/tot_sum.loc[solid], 
               c=cmap['green'], s=10, marker='o', alpha=1)
    ax.scatter([dates_axis.index(x) for x in faded], datedf.loc[faded].refpos/tot_sum.loc[faded], 
               c=cmap['green'], s=10, marker='o', alpha=.35)

    idx_eu = metadf.loc[idx_ok][metadf.loc[idx_ok].region=='Europe'].index
    datedf = pd.concat([metadf.loc[idx_eu].date, mut_freqdf.loc(axis=0)[idx_eu, :].sum(level='sample')],axis=1)
    datedf = datedf.set_index('date').sum(level='date')
    add_dates = set(dates_axis).difference(datedf.index)
    add_dates = pd.DataFrame([np.nan]*len(add_dates), index=add_dates, columns=['refpos'])
    datedf = datedf.append(add_dates)
    datedf.sort_index(inplace=True)
    tot_sum = metadf.loc[idx_eu].groupby('date').count().acc.append(add_dates.refpos)
    solid = [date for (date, val) in zip(dates_axis, tot_sum.loc[dates_axis].values) if val>=5]
    faded = list(set(dates_axis).difference(solid))
    ax.scatter([dates_axis.index(x) for x in solid], datedf.loc[solid].refpos/tot_sum.loc[solid], 
               c=cmap['orange'], s=10, marker='o', alpha=1)
    ax.scatter([dates_axis.index(x) for x in faded], datedf.loc[faded].refpos/tot_sum.loc[faded], 
               c=cmap['orange'], s=10, marker='o', alpha=.35)
    
    datedf = pd.concat([metadf.loc[idx_ok].date, mut_freqdf.loc(axis=0)[idx_ok, :].sum(level='sample')],axis=1)
    datedf = datedf.set_index('date').sum(level='date')
    add_dates = set(dates_axis).difference(datedf.index)
    add_dates = pd.DataFrame([np.nan]*len(add_dates), index=add_dates, columns=['refpos'])
    datedf = datedf.append(add_dates)
    datedf.sort_index(inplace=True)
    tot_sum = metadf.loc[idx_ok].groupby('date').count().acc.append(add_dates.refpos)
    solid = [date for (date, val) in zip(dates_axis, tot_sum.loc[dates_axis].values) if val>=5]
    faded = list(set(dates_axis).difference(solid))
    ax.scatter([dates_axis.index(x) for x in solid], datedf.loc[solid].refpos/tot_sum.loc[solid], 
               c=cmap['blue'], s=10, marker='o', alpha=1)
    ax.scatter([dates_axis.index(x) for x in faded], datedf.loc[faded].refpos/tot_sum.loc[faded], 
               c=cmap['blue'], s=10, marker='o', alpha=.35)
    
    ax.set_xlim(-1, len(dates_axis)+.5)
    ax.set_ylabel('Number of mutations per sample')
    ax.grid('on', alpha=.5)
    ax.xaxis.set_ticks(range(0, len(dates_axis)+1, 7))
    ax.xaxis.set_ticklabels(pd.date_range(start=dates_axis[::7][0], freq='W', periods=len(ax.get_xticks())).format('%D-%M-%Y')[1:])
    ax.xaxis.set_tick_params(labelrotation=90)
    ax.set_ylim(0, 14.65)
    
    leg = [Line2D([0],[0], color='w', label='Netherlands', marker='o', markerfacecolor=cmap['green'], markersize=10),
           Line2D([0],[0], color='w', label='Europe', marker='o', markerfacecolor=cmap['orange'], markersize=10),
           Line2D([0],[0], color='w', label='Global', marker='o', markerfacecolor=cmap['blue'], markersize=10)]
    plt.legend(handles=leg, loc='upper left', fontsize=10)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig3.png'), dpi=600, format='png')
    plt.savefig(os.path.join(output_dir, 'fig3.pdf'), dpi=600, format='pdf')


