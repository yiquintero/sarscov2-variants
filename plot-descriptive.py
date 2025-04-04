#!/usr/bin/env python3
"""
This script generates the figures .... from our publication
It needs to location of SARS-CoV-2 metadata file and the output directory for images
"""
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from datetime import datetime

def main():
    parser = argparse.ArgumentParser(description="Plot descriptive statistics of the total SARS-CoV-2 population and the individual clades)
    parser.add_argument("metadata", help="SARS-CoV-2 metadata file")
    parser.add_argument('-o', '--output-dir', type=Path, required=True,
                        help="Output directory to save the images.")

    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    metadf = pd.read_csv(metadata, index_col=0, header=0, sep='\t')

    sns.set(style="whitegrid", palette="muted", color_codes=True)
    plt.figure(figsize=(7, 6))
    sns.barplot(x=sorted(metadf.region.unique()), y=metadf.groupby('region').count().date, 
                color='b')
    plt.xlabel('Region')
    plt.ylabel('Number of genomes')
    plt.grid('on')
    plt.savefig('img/desc-totalregion_all.svg', dpi=600, format='svg')

if __name__ == "__main__":
    main()
