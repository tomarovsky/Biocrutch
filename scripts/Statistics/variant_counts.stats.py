#!/usr/bin/env python3
__author__ = 'tomarovsky'
from argparse import ArgumentParser
import pandas as pd


def main():
    for file in args.input:
        df = pd.read_csv(file, sep='\t', header=0)
        sample = df.columns[2]

        if args.exclude:
            df = df[df['CHROM'] != args.exclude]

        df[sample] = df[sample] / args.window_size * args.multiplicator

        mean = df[sample].mean()
        median = df[sample].median()

        print(f'{sample}\t{mean:.2f}\t{median:.2f}')


if __name__ == '__main__':
    parser = ArgumentParser(description='script to calculate mean and median based .variant_counts.tsv (SNPs/1kbp)')
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str, nargs='+', help='MACE BED file')
    group_required.add_argument('-w', '--window_size', type=int, help='Window size', default=1000000)
    group_required.add_argument('-m', '--multiplicator', type=int, help='Multiplicator', default=1000)
    group_additional = parser.add_argument_group('Additional options')
    group_additional.add_argument('-e', '--exclude', type=str, default=False, help='Exclude chr')
    args = parser.parse_args()
    main()
