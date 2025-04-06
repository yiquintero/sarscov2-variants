#!/usr/bin/env python3
"""
This script generates Figure 7 from our publication
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
    parser = argparse.ArgumentParser(description="Plot mutation frequencies of SARS-CoV-2 genomes in the Netherlands "
                                     "and the total number of genomes submitted to the databases")
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
    countries = ['UK', 'USA', 'Australia', 'Canada', 'India', 'Spain', 'China', 'Netherlands']
    nlargest = 15
    maxp = len(countries)

    warnings.simplefilter(action='ignore') 

    print("Plotting Figure 7")
    fig, ax = plt.subplots(2, 1, figsize=(12, 6))

    # We plot only the samples that have the full date available
    metadf.loc[:, 'acc'] = metadf.index # store the genome accession IDs for future use
    idx_all = metadf[metadf.date.apply(lambda x: len(x)==10)].index.intersection(set(mutdf.index))
    all_variants = [sorted(mutdf.loc[metadf.loc[idx_all][metadf.loc[idx_all].country==c].index].groupby('varname').count().refpos.nlargest(nlargest).index)
                    for c in countries]
    uniq_variants, variant_counts = np.unique(np.array(all_variants), return_counts=True)
    variant_count_dict = dict(zip(uniq_variants, variant_counts))

    freqdf = mutdf.groupby(['sample','varname']).count().refpos
    dates_axis = sorted(metadf.loc[idx_all].date.unique())
    dates_axis = pd.date_range(start=dates_axis[0], end=dates_axis[-1], freq='D').format('%D-%M-%Y')[1:]

    country = 'Netherlands'
    idx_nl = metadf.loc[idx_all][metadf.loc[idx_all].country==country].index
    toplot = mutdf.loc[idx_nl].groupby('varname').count().refpos.nlargest(nlargest)
    for i, p in enumerate(toplot.index):
        datedf = pd.concat([metadf.loc[idx_nl].date, freqdf.loc(axis=0)[idx_nl, p].sum(level='sample')],axis=1)
        datedf = datedf.set_index('date').sum(level='date')
        add_dates = set(dates_axis).difference(datedf.index)
        add_dates = pd.DataFrame([np.nan]*len(add_dates), index=add_dates, columns=['refpos'])
        datedf = datedf.append(add_dates)
        datedf.sort_index(inplace=True)
        tot_sum = metadf.loc[idx_nl].groupby('date').count().acc.append(add_dates.refpos)
        color = '#8c8c8c'
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
        ax[0].plot(dates_axis, ma_perc, color=color)
    ax[0].set_ylabel('Frequency')
    ax[0].set_ylim(0, 1)
    ax[0].set_title(country)
    ax[0].grid('on', alpha=.5)
    ax[0].xaxis.set_ticks(dates_axis[::7])
    ax[0].set_xlim(56, 163)
    ax[0].xaxis.set_ticklabels([])

    variant2color = {'Shared': '#4978d0', 'Unique': '#d66060'}
    leg = [Line2D([0],[0], color='w', label=k, marker='s', markerfacecolor=v, markersize=5) 
           for k,v in variant2color.items()]
    ax[0].legend(handles=leg, loc='center left', fontsize=8)
    tick_labels = pd.date_range(start=dates_axis[::7][0], freq='W', periods=len(ax[0].get_xticks())).format('%D-%M-%Y')[1:]

    w = 7
    minp = 5
    sub_metadf = metadf[metadf.date.apply(lambda x: len(x)==10)][['country','date-fixed', 'recName','region','clade',
                                                                  'lineage','age','subloc','timestmp','numDupes',
                                                                  'gender','gisaid-clade','comp-clade','LS-clade']]
    country2color = {'Australia': '#4978d0', 'Netherlands': '#ee854b', 'UK': '#6acc65', 
                    'USA': '#956bb3', 'Canada': '#dd7ec0', 'India': '#d5bb67'}
    countries = sorted(['Australia', 'Netherlands', 'UK', 'USA', 'Canada', 'India'])
    datedf = sub_metadf[sub_metadf.country.apply(lambda x: x in countries)].groupby(['date-fixed','country']).count().copy()
    datedf.sort_index(level='date-fixed', inplace=True)
    datedf.sort_index(level='country', inplace=True)

    datedf = sub_metadf[sub_metadf.country.apply(lambda x: x in countries)].set_index('date-fixed')
    datedf.sort_index(axis=0, inplace=True)
    toplot = datedf[datedf.country==country].groupby(by='date-fixed').count().recName
    add_date = set(dates_axis).difference(set(toplot.index))
    toplot = toplot.append(pd.Series([0]*len(add_date), index=add_date))
    toplot.sort_index(inplace=True)
    ma = pd.DataFrame(toplot.loc[dates_axis]).rolling(w, min_periods=minp).mean()
    ma[ma==0] = np.nan
    ax[1].plot(dates_axis, ma.loc[dates_axis], color='k')
    ax[1].set_ylabel('# of genomes')
    ax[1].grid('on', alpha=.5)
    ax[1].xaxis.set_ticks(dates_axis[::7])
    ax[1].xaxis.set_ticklabels(dates_axis[::7])
    ax[1].xaxis.set_tick_params(labelrotation=90)
    ax[1].set_xlabel('Date')
    ax[1].set_xlim(56,163)
    ax[1].set_ylim(0,30)
    ax[0].set_title('')
    fig.set_size_inches(9.43,4.4)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig7.png'), dpi=600, format='png')
    plt.savefig(os.path.join(output_dir, 'fig7.pdf'), dpi=600, format='pdf')

if __name__ == "__main__":
    main()    