#!/usr/bin/env python3
__author__ = 'tomarovsky'
'''
Script for calculating median, average, maximum and minimum coverage.
a. calculate stats for whole genome
b. calculate stats for each scaffold
c. calculate stats in 100 kbp and 1 Mbp stacking windows
Header-less tab-separated input file with 3 columns: scaffold_id, position(1-based), coverage.
'''
from Biocrutch.Statistics.coverage_statistics.get_coverage_statistics import GetCoverageStatistics
from Biocrutch.Statistics.coverage_statistics.coverage_metrics import CoveragesMetrics
from Biocrutch.Routines.routine_functions import metaopen
from sys import stdin
import pandas as pd
import argparse


def main():
    metrics = GetCoverageStatistics(args.input, args.output)

    if args.all_stats:
        metrics.get_all_statistics(args.frame_size)
    if args.whole_genome_stats:
        metrics.get_whole_genome_stats()
    if args.scaffolds_stats:
        metrics.get_scaffolds_stats()
    if args.nonoverlapping_windows_stats:
        metrics.get_nonoverlapping_windows_stats(args.frame_size)
    if args.overlapping_windows_stats: # in developing
        metrics.get_overlapping_windows_stats(args.frame_size)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="script for calculating median, average, maximum and minimum coverage in file.tab. Report to files.csv")

    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=lambda s: metaopen(s, "rt"),
                                  help="input file.bam.gz (don`t use for STDIN)", default=stdin)

    group_additional = parser.add_argument_group('Additional options')
    group_additional.add_argument('-f', '--frame-size', type=int,
                                help="calculate stats in 100 kbp and 1 Mbp stacking windows", default=1000000)
    group_additional.add_argument('-o', '--output', metavar='PATH', type=str, default=False,
                                help='output file prefix without frame size')
    group_additional.add_argument('--all-stats', action="store_true", default=False,
                                help="to calculate statistics for overlapping and non-overlapping windows, scaffolds and whole genome")
    group_additional.add_argument('--whole-genome-stats', action="store_true", default=False,
                                help="to calculate statistics for whole genome only")
    group_additional.add_argument('--scaffolds-stats', action="store_true", default=False,
                                help="to calculate statistics for scaffolds only")
    group_additional.add_argument('--nonoverlapping-windows-stats', action="store_true", default=False,
                                help="to calculate statistics for non-overlapping windows only")
    group_additional.add_argument('--overlapping-windows-stats', action="store_true", default=False,
                                help="to calculate statistics for overlapping windows only")

    args = parser.parse_args()
    main()
