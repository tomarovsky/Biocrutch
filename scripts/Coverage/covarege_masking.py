#!/usr/bin/env python3
__author__ = 'tomarovsky'
'''
Get coverage mask. For pipeline:
zcat SAMPLE.per-base.bed.gz | awk '{if( ($4 > 2.5 * whole_genome_median) || ($4 < 0.5 * whole_genome_median) print $0)}' | covarege_masking.py
zcat SAMPLE.per-base.bed.gz | awk '{if( ($4 > 2.5 * $(cat *_whole_genome_stats.csv | sed -n 2p | awk '{print $2}')) || ($4 < 0.5 * $(cat *_whole_genome_stats.csv | sed -n 2p | awk '{print $2}')) print $0)}' | covarege_masking.py
'''
from Biocrutch.Routines.routine_functions import metaopen
import argparse
from sys import stdin


def main():
    barcode_fd = metaopen(args.prefix + "_barcodes.txt.gz", 'wt')
    forward_fd = metaopen(args.prefix + "_ema-bin-all_1.fastq.gz", 'wt')
    reverse_fd = metaopen(args.prefix + "_ema-bin-all_2.fastq.gz", 'wt')
    if args.input == '-':
        data = stdin
    elif args.compressed: #for file.gz
        data = metaopen(args.input, "rt", args.buffering)
    else:
        data = metaopen(args.input, "r", args.buffering)
    for line in data:
        line = line.split()
        barcode_fd.write(line[0] + '\n')
        forward_fd.write(line[1] + '\n' + line[2] + '\n' + '+\n' + line[3] + '\n')
        reverse_fd.write(line[1] + '\n' + line[4] + '\n' + '+\n' + line[5] + '\n')
    barcode_fd.close()
    forward_fd.close()
    reverse_fd.close()
    if args.input != '-':
        data.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="get the coverage mask through the pipeline")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=lambda s: metaopen(s, "rt"),
                                help="input coverage_statistics_output.csv (don`t use for STDIN)", default=stdin)
    group_required.add_argument('-o', '--output', type=str,
                                help="output file prefix")
    args = parser.parse_args()
    main()