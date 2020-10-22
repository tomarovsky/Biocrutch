#!/usr/bin/env python3
__author__ = 'tomarovsky'
'''
Script for determining the coordinates of the pseudoautosomal region.
Output of coordinates to BED file.
'''
from Biocrutch.Routines.routine_functions import metaopen
from sys import stdin
import pandas as pd
import argparse

# function to convert output to BED format. Accepts most likely a dictionary


def main():
    deviation_percent = 25 # the percent of deviation
    deviation = args.whole_genome_value / 100 * deviation_percent
    minimum_coverage = args.whole_genome_value - deviation # 25.5
    maximum_coverage = args.whole_genome_value + deviation # 42.5
    count_repeat_frames = 0

    for line in args.input:
        line = line.rstrip().split("\t")
        coverage_value = int(line[args.coverage_column_name])
        if coverage_value > minimum_coverage and coverage_value < maximum_coverage:
            count_repeat_frames += 1
        else:
            count_repeat_frames = 0
        print(count_repeat_frames)
        




        


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script for determining the coordinates of the pseudoautosomal region. Output of coordinates to BED file.")

    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=lambda s: metaopen(s, "rt"),
                                  help="input file.bam.gz (don`t use for STDIN)", default=stdin)

    group_additional = parser.add_argument_group('Additional options')
    group_additional.add_argument('-o', '--output', metavar='PATH', type=str, default=False,
                                  help='output file prefix')
    group_additional.add_argument('-f', '--frame-size', type=int,
                                  help="calculate stats in 100 kbp and 1 Mbp stacking windows", default=100000)
    group_additional.add_argument('--coverage_column_name', type=int,
                                  help="Number of column in coverage file with mean/median coverage per window", default=3)
    group_additional.add_argument('-s', '--scaffold-name', type=str,
                                  help="Name of column in coverage file with scaffold name", default="scaffold")
    group_additional.add_argument('-m', '--whole_genome_value', type=int,
                                  help="whole genome median/mean value", default=34)

    args = parser.parse_args()
    main()