#!/usr/bin/env python3
__author__ = 'tomarovsky'
"""
Get coverage mask. You can use pipeline:
zcat SAMPLE.per-base.bed.gz | awk '{if( ($4 > 2.5 * whole_genome_median) || ($4 < 0.5 * whole_genome_median) print $0)}' | coverage_masking.py
WHOLE_MEDIAN=$(cat *_whole_genome_stats.csv | sed -n 2p | awk '{print $2}');
zcat SAMPLE.per-base.bed.gz | awk '{if( ($4 > 2.5 * $WHOLE_MEDIAN) || $4 < 0.5 * $WHOLE_MEDIAN) print $0)}' | coverage_masking.py
"""
from Biocrutch.Routines.routine_functions import metaopen, metaoutput
import argparse
from sys import stdin


def main():
    outfile = metaopen(metaoutput(args.output, '.mask.bed.gz'), 'wt')
    for line in args.input:
        scaffold_name, start, stop, coverage_value = args.input.readline().strip().split()
        coverage_value = int(coverage_value)
        if (coverage_value > (2.5 * args.whole_median)) or coverage_value < (0.5 * args.whole_median):
            break

    for i in args.input:
        line = i.strip().split()
        line[3] = int(line[3])
        if (line[3] > (2.5 * args.whole_median)) or (line[3] < (0.5 * args.whole_median)):
            if line[0] == scaffold_name:
                if stop == line[1]:
                    stop = line[2]
                    continue
                else:
                    outfile.write('\t'.join([scaffold_name, start, stop]) + '\n')
                    start = line[1]
                    stop = line[2]
            else:
                outfile.write('\t'.join([scaffold_name, start, stop]) + '\n')
                scaffold_name, start, stop = line[:3]
    outfile.write('\t'.join([scaffold_name, start, stop]) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='get the coverage mask through the pipeline')
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=lambda s: metaopen(s, 'rt'), help='input coverage_statistics_output.csv)', default=stdin)
    group_required.add_argument('-w', '--whole-median', type=float, help='whole median value')
    group_required.add_argument('-o', '--output', type=str, help='output file prefix')
    args = parser.parse_args()
    main()
