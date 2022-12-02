#!/usr/bin/env python3
import matplotlib.pyplot as plt
from venn import venn
plt.ioff()
import pandas as pd
import argparse



def df_csv(path):
    with open(path) as file:
        comment_lines = 0
        for line in file:
            if line.startswith('#'):
                comment_lines += 1
                names = line.replace('\n', '').split('\t')
            else:
                comment_lines -= 1
                break
    df = pd.read_csv(path, header=0, names=names, skiprows=comment_lines, sep='\t')
    return df


def main():
    statuses = ['Complete', 'Duplicated', 'Fragmented', 'Missing']
    df_list = [df_csv(file) for file in args.full_tables]
    if not args.species:
        species = [file.split('_')[-1][:-4] for file in args.full_tables] # if not args.species else args.species
    else:
        species = args.species
    for status in statuses:
        fig, ax = plt.subplots(figsize=(16, 8), nrows=1, ncols=2)
        df_status_list = [set(list(df.loc[df['Status'] == status]['# Busco id'])) for df in df_list]
        venn_dict = dict(zip(species, df_status_list))
        venn(venn_dict, ax=ax[0], fontsize=11)
        venn(venn_dict, ax=ax[1],fontsize=11, fmt="{percentage:.1f}%")
        ax[0].set_title(status)
        ax[1].set_title(status)
        plt.savefig(f'{status}.venn.png', dpi = 500)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Venn diagrams")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-f', '--full-tables', type=str, nargs='+',
                                help="full_table_SPECIESNAME.tsv files")
    group_required.add_argument('-s', '--species', default=False, type=str,
                                nargs='+', help="full_table_SPECIESNAME.tsv files")
    group_required.add_argument('-o', '--outdir', type=str, help="output directory name")
    args = parser.parse_args()
    main()