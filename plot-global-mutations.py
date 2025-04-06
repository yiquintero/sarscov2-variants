#!/usr/bin/env python3
"""
This script generates the figures 5 and 6 from our publication
It needs to location of SARS-CoV-2 metadata file, mutation file (obtained from coronapp)
and the output directory for images
"""
import os
import argparse
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Plot the total frequency (\% of genomes) of the top 15 mutations in the"
                                                 " most-sampled countries in our dataset and how the freuqencies of these"
                                                 " mutations change over time")
    parser.add_argument("metadata", help="SARS-CoV-2 metadata file")
    parser.add_argument('mutation_file', help="SARS-CoV-2 mutations from coronapp")
    parser.add_argument('-o', '--output-dir', type=Path, required=True,
                        help="Output directory to save the images.")
    args = parser.parse_args()
    output_dir = args.output_dir
    metadata_file = args.metadata
    mutation_file = args.mutation_file
    os.makedirs(output_dir, exist_ok=True)

    # Load metadata file describing each SARS-CoV-2 genome
    metadf = pd.read_csv(metadata_file, index_col=0, header=0, sep='\t')
    mutdf = pd.read_csv(mutation_file, sep='\t', index_col=0, header=0)

    # Define which countries to plot the variants for
    countries = ['Netherlands', 'UK', 'USA', 'Australia', 'Canada', 'India', 'Spain', 'China']
    variant2color = {'Shared': '#4978d0', 'Non-unique non-common': 'gray', 'Unique': '#d66060'}
    nlargest = 15
    # We plot only the samples that have the full date available
    idx_all = metadf[metadf.date.apply(lambda x: len(x)==10)].index.intersection(set(mutdf.index))
    all_variants = [sorted(mutdf.loc[metadf.loc[idx_all][metadf.loc[idx_all].country==c].index].groupby('varname').count().refpos.nlargest(nlargest).index)
                    for c in countries]
    uniq_variants, variant_counts = np.unique(np.array(all_variants), return_counts=True)
    variant_count_dict = dict(zip(uniq_variants, variant_counts))
    
    # Figure 5
    warnings.simplefilter(action='ignore')
    print("Plotting Figure 5")
    fig, ax = plt.subplots(2, 4, figsize=(12, 8))
    for xy,country in zip(sum(ax.tolist(),[]),countries):
        idx_country = metadf.loc[idx_all][metadf.loc[idx_all].country==country].index
        toplot = mutdf.loc[idx_country].groupby('varname').count().refpos.nlargest(15)/len(idx_country)
        color = ['gray']*15
        for k,v in enumerate(toplot.index):
            if variant_count_dict[v]==1:
                color[k] = '#d66060'
            elif variant_count_dict[v] == len(countries):
                color[k] = '#4978d0'
        xy.bar(toplot.index, toplot.values, color=color)
        xy.xaxis.set_ticklabels(toplot.index)
        xy.xaxis.set_tick_params(rotation=90)
        xy.set_title(country)
        xy.grid('on', alpha=.5)
        xy.set_xlim(-1, 15)
        xy.set_ylabel('% of genomes')
    
    ax[1, 3].set_ylim(0, 0.8)  
    leg = [Line2D([0],[0], color='w', label=k, marker='s', markerfacecolor=v, markersize=10)
           for k,v in variant2color.items()]
    plt.legend(handles=leg, loc='upper right', fontsize=10)
    for xy in sum(ax.tolist(),[])[:-1]: 
        xy.set_ylim(0, 0.86)
    for i in [1, 2, 3]:
        ax[0, i].set_ylabel('')
        ax[1, i].set_ylabel('')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig5.png'), dpi=600, format='png')
    plt.savefig(os.path.join(output_dir, 'fig5.pdf'), dpi=600, format='pdf')

    # Figure 6
    metadf.loc[:, 'acc'] = metadf.index # store the genome accession IDs for future use
    nlargest = 15
    # Making sure the dates match
    mutdf.loc[:, 'date'] = metadf.loc[mutdf.index.intersection(metadf.index), 'date']
    countries = ['UK', 'USA', 'Australia', 'Canada', 'India', 'Spain']
    maxp = len(countries)
    idx_all = metadf[metadf.date.apply(lambda x: len(x)==10)].index.intersection(set(mutdf.index))
    all_variants = [sorted(mutdf.loc[metadf.loc[idx_all][metadf.loc[idx_all].country==c].index].groupby('varname').count().refpos.nlargest(nlargest).index)
                    for c in countries]
    uniq_variants, variant_counts = np.unique(np.array(all_variants), return_counts=True)
    variant_count_dict = dict(zip(uniq_variants, variant_counts))
    
    mut_freqdf = mutdf.groupby(['sample','varname']).count().refpos
    dates_axis = sorted(metadf.loc[idx_all].date.unique())
    dates_axis = pd.date_range(start=dates_axis[0], end=dates_axis[-1], freq='D').format('%D-%M-%Y')[1:]
    
    plt.style.use('seaborn-muted')
    fig, ax = plt.subplots(3,2, figsize=(12,6))
    for xy,country in zip(sum(ax.tolist(),[]), countries):
        idx_country = metadf.loc[idx_all][metadf.loc[idx_all].country==country].index
        toplot = mutdf.loc[idx_country].groupby('varname').count().refpos.nlargest(nlargest)
        for i,p in enumerate(toplot.index):
            datedf = pd.concat([metadf.loc[idx_country].date, mut_freqdf.loc(axis=0)[idx_country, p].sum(level='sample')],axis=1)
            datedf = datedf.set_index('date').sum(level='date')
            add_dates = set(dates_axis).difference(datedf.index)
            add_dates = pd.DataFrame([np.nan]*len(add_dates), index=add_dates, columns=['refpos'])
            datedf = datedf.append(add_dates)
            datedf.sort_index(inplace=True)
            tot_sum = metadf.loc[idx_country].groupby('date').count().acc.append(add_dates.refpos)
            color = '#8c8c8c'
            # Set different minimum and/or maximum period limit to make the plot prettier
            if variant_count_dict[p]==1:
                color = '#d66060'
            elif variant_count_dict[p] == maxp:
                color = '#4978d0'
            if country == 'USA':
                minp = 7
            elif country == 'Canada':
                minp = 0
            else:
                minp = 5
            ma_perc = (datedf.loc[dates_axis].refpos/tot_sum.loc[dates_axis]).rolling(7, min_periods=minp).mean()
            xy.plot(dates_axis, ma_perc, color=color)
        xy.set_xlim(56, 163)
        xy.set_ylabel('Frequency')
        xy.set_ylim(0, 1)
        xy.set_title(country)
        xy.grid('on', alpha=.5)
        xy.xaxis.set_ticks(dates_axis[::7])
        xy.xaxis.set_ticklabels([])
    
    leg = [Line2D([0],[0], color='w', label=k, marker='s', markerfacecolor=v, markersize=10)
           for k,v in variant2color.items()]
    tick_labels = pd.date_range(start=dates_axis[::7][0], freq='W', periods=len(xy.get_xticks())).format('%D-%M-%Y')[1:]
    ax[2, 1].legend(handles=leg, loc='center right', fontsize=10)
    ax[2, 1].xaxis.set_ticklabels(tick_labels)
    ax[2, 1].xaxis.set_tick_params(labelrotation=90)
    ax[2, 0].xaxis.set_ticklabels(tick_labels)
    ax[2, 0].xaxis.set_tick_params(labelrotation=90)
    ax[2, 0].set_xlabel('Date')
    ax[2, 1].set_xlabel('Date')
    for i in [0, 1, 2]: 
        ax[i,1].set_ylabel('')
    ax[0, 0].set_xlim(56, 163)
    ax[0, 1].set_xlim(56, 163)
    ax[1, 1].set_xlim(56, 163)
    ax[1, 0].set_xlim(56, 163)
    ax[2, 0].set_xlim(56, 163)
    ax[2, 1].set_xlim(56, 163)
    fig.set_size_inches(9.43, 6.6)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig6.png'), dpi=600, format='png')
    plt.savefig(os.path.join(output_dir, 'fig6.pdf'), dpi=600, format='pdf')
