#!/usr/bin/env python3
__author__ = 'tomarovsky'
from argparse import ArgumentParser, FileType
from Biocrutch.Parsers.fasta_opener import Fasta_opener
from collections import defaultdict
from Biocrutch.Statistics.quast_core.quast_stats import Quast_core
import Biocrutch.Statistics.quast_core.constants as c
import pandas as pd
from os import path


def main():
    df = pd.DataFrame()
    data_dict = defaultdict(list)
    for file in args.input:
        data = Fasta_opener(file)
        seq_dict = data.parse_sequences(args.buffering)
        lengths = data.lengths_to_frame()
        metrics = Quast_core(seq_dict, lengths)

        for min_contig in args.min_contig:
            data_dict['Contigs(>={})'.format(str(min_contig))].append(metrics.contig_count(min_contig))

        for min_contig in args.min_contig:
            data_dict['Totallength(>={})'.format(str(min_contig))].append(metrics.total_length(min_contig))

        for min_contig in args.min_contig:
            data_dict['GC(%)(>={})'.format(str(min_contig))].append(metrics.gc_content(min_contig))

        for min_contig in args.min_contig:
            data_dict['N_amount(>={})'.format(str(min_contig))].append(metrics.n_amount(min_contig))

        data_dict['Largestcontig'].append(metrics.largest_contig_lengh())

        for i in args.nl_statistics:
            for min_contig in args.min_contig:
                n_l_stat = metrics.n_l_statistics(i, min_contig)
                data_dict['N{0}stats(>={1})'.format(str(i), str(min_contig))].append(n_l_stat[0])
                data_dict['L{0}stats(>={1})'.format(str(i), str(min_contig))].append(n_l_stat[1])

    df = pd.DataFrame(data_dict, index=args.input)

    if args.print:  # print to terminal
        print(df.T.astype(str))
    if args.output:  # create a report.csv
        df.to_csv(args.output, encoding='utf-8')

    return df


if __name__ == '__main__':
    parser = ArgumentParser(description=c.DESCRIPTION)

    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str, nargs='+', required=True, help=c.HELP_INPUT)

    group_additional = parser.add_argument_group('Additional options')
    group_additional.add_argument('-b', '--buffering', metavar='INT', type=int, default=None, help=c.HELP_BUFFERING)
    group_additional.add_argument('-m', '--min-contig', metavar='INT', type=int, nargs='+', default=[0, 150], help=c.HELP_MIN_CONTIG)
    group_additional.add_argument('-n', '--nl-statistics', metavar='INT', type=int, nargs='+', default=[50, 75], help=c.HELP_NL_STATISTICS)
    group_additional.add_argument('-o', '--output', type=str, help=c.HELP_OUTPUT)
    group_additional.add_argument('-p', '--print', type=lambda x: (str(x).lower() in ['true', 'yes', '1']), default=True, help=c.HELP_PRINT)
    args = parser.parse_args()

    main()
