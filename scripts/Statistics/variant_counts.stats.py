#!/usr/bin/env python3
__author__ = 'tomarovsky'
from argparse import ArgumentParser
import pandas as pd


def main():
    for file in args.input:
        df = pd.read_csv(args.input, sep='\t', header=0)
        sample = df.columns[2]

        # Calculate mean, median and mode
        mean = df[sample].mean()
        median = df[sample].median()
        mode = df[sample].mode()[0]

        print(f"{sample}\t{mean}\t{median}\t{mode}")


if __name__ == "__main__":
    parser = ArgumentParser(description="script to calculate mean, median and mode based .variant_counts.tsv")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str, nargs='+', help="MACE BED file")
    # group_additional = parser.add_argument_group('Additional options')
    # group_additional.add_argument('-o', '--output', type=str, default=False, help="outfile name")
    # group_additional.add_argument('-h', '--header', type=bool, default=False, help="write header to outfile")
    args = parser.parse_args()
    main()

