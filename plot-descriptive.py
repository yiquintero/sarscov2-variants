#!/usr/bin/env python3
"""
This script generates the figures 1A, 1B and 2 from our publication
It needs to location of SARS-CoV-2 metadata file and the output directory for images
"""
import os
import argparse
import warnings
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Plot descriptive statistics of "
                                     "the total SARS-CoV-2 population and the individual clades")
    parser.add_argument("metadata", help="SARS-CoV-2 metadata file")
    parser.add_argument('-o', '--output-dir', type=Path, required=True,
                        help="Output directory to save the images.")

    args = parser.parse_args()
    metadata_file = args.metadata
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    # Load metadata file describing each SARS-CoV-2 genome
    metadf = pd.read_csv(metadata_file, index_col=0, header=0, sep='\t')

    # Figure 1
    print("Plotting Figure 1")
    warnings.simplefilter(action='ignore')
    sns.set_theme(style="whitegrid", palette="muted", color_codes=True)
    fig = plt.figure(figsize=(12, 6))
    # Figure 1A
    ax1 = fig.add_subplot(1, 2, 1)
    sns.barplot(x=sorted(metadf.region.unique()), y=metadf.groupby('region').count().date, 
                color='b', ax=ax1)
    ax1.set_xlabel('Region')
    ax1.set_ylabel('Number of genomes')
    ax1.grid('on')

    # Figure 1B
    ax2 = fig.add_subplot(1, 2, 2)
    new_metadf = metadf.copy()
    # Cleaning the "date" field in metadata file plot correctly
    new_metadf.loc[:, 'date'] = new_metadf['date-fixed'].apply(lambda x: '-'.join(x.split('-')[:2]))
    new_metadf.date.replace({'2020-5': '2020-05'}, inplace=True)
    new_metadf = new_metadf.groupby(['date', 'region']).count()
    new_metadf = new_metadf.sort_index(level=0)
    dates_axis = new_metadf.index.get_level_values('date').unique()

    for region in metadf.region.unique():
        # Padding the data from each region to have matching x axes
        add_dates = dates_axis.difference(new_metadf.loc(axis=0)[:, region].index.get_level_values('date'))
        for d in add_dates: 
            new_metadf.loc(axis=0)[d, region] = 0

    plt.style.use('seaborn-muted')
    for region in sorted(metadf.region.unique()):
        a = new_metadf.loc(axis=0)[:, region].apply(np.cumsum)
        ax2.plot(dates_axis, a.country, markersize=4, marker='o')
    ax2.set_xlabel("Date", fontsize= 12)
    ax2.set_ylabel("Number of genomes", fontsize= 12)
    ax2.set_xlim([min(dates_axis), max(dates_axis)])
    ax2.set_ylim(ax1.get_ylim())
    ax2.legend(sorted(metadf.region.unique()), loc='upper left')
    
    ax1.text(-1, 17000, 'A',fontsize=15)
    ax2.text(-1, 17000, 'B',fontsize=15)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig1.png"), dpi=600, format='png')
    plt.savefig(os.path.join(output_dir, "fig1.pdf"), dpi=600, format='pdf')

    # Figure 2
    print("Plotting Figure 2")
    plt.style.use('seaborn-muted')
    clean_metadf = metadf[metadf.date.apply(lambda x: len(x)==10)]
    clean_metadf = clean_metadf[['country','date-fixed', 'recName','region','clade','lineage',
                                 'age','subloc','timestmp','numDupes','gender',
                                 'gisaid-clade','comp-clade','LS-clade']]
    # Set moving window size to 7 days = 1 week 
    w = 7 
    # Cleaning the dates axis to plot prettier
    dates_axis = sorted(clean_metadf['date-fixed'].unique())
    dates_axis_ticklabel = dates_axis[w-1::7]
    dates_axis_ticklabel[0] = '2020-01-08'
    clade2color = {'L': '#4978d0', 'S': '#ee854b', 'G': '#6acc65', 'GR': '#956bb3',
                  'GH': '#dd7ec0', 'V': '#d5bb67'}
    countries = ['Netherlands', 'UK', 'Australia', 'Singapore', 'China']
    fig, ax = plt.subplots(5, 1, figsize=(12, 32))
    
    for i, country in enumerate(countries):
        clean_metadf_fig = clean_metadf[clean_metadf.country==country].groupby(['date-fixed','gisaid-clade']).count()
        clean_metadf_fig.sort_index(level='date-fixed', inplace=True)
        clades = clean_metadf_fig.groupby('gisaid-clade').sum().index
        for clade in clades:
            add_dates = set(dates_axis).difference(clean_metadf_fig.loc(axis=0)[:, clade].index.get_level_values('date-fixed'))
            for d in add_dates: 
                clean_metadf_fig.loc(axis=0)[d, clade] = 0
        clean_metadf_fig = clean_metadf_fig.sort_index(level='date-fixed')
        total_sum = clean_metadf_fig.sum(level='date-fixed').country.values.reshape(1, len(dates_axis))
        clade_sum = clean_metadf_fig.country.values.reshape(len(dates_axis), len(clades)).T
        ma_clade_sum = pd.DataFrame((clade_sum/total_sum*100).T).rolling(w, min_periods=0).mean().values.T
        colors = [clade2color[clade] for clade in clades]
        ax[i].grid('y', alpha=.5)
        ax[i].stackplot(dates_axis[w-1:], ma_clade_sum[:, w-1:], colors=colors, labels=clades)
        ax[i].set_ylabel("% of genomes", fontsize=12)
        ax[i].xaxis.set_ticks([])
        ax[i].xaxis.set_tick_params(grid_alpha=0)
        ax[i].set_xlim([dates_axis[w-1::7][0], dates_axis[w-1::7][-1]])
        ax[i].set_ylim([0, 100])
        ax[i].set_title(country)

    ax[i].xaxis.set_ticks(dates_axis[w-1::7])
    ax[i].xaxis.set_ticklabels(dates_axis_ticklabel)
    ax[i].xaxis.set_tick_params(labelrotation=90, grid_alpha=0)
    ax[i].legend(clades, loc='upper right')
    plt.tight_layout()

    plt.savefig(os.path.join(output_dir, "fig2.png"), dpi=600, format='png')
    plt.savefig(os.path.join(output_dir, "fig2.pdf"), dpi=600, format='pdf')
    
if __name__ == "__main__":
    main()
