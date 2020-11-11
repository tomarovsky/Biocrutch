#!/usr/bin/env python3
__author__ = 'tomarovsky'
'''
Script for determining the coordinates of the pseudoautosomal region.
Output of coordinates to BED file.
'''
from Biocrutch.Statistics.coverage_statistics.coverage_metrics import CoveragesMetrics
from Biocrutch.Statistics.pseudoautosomal_region.coordinates import Coordinates
from Biocrutch.Routines.routine_functions import metaopen
from collections import Counter
from sys import stdin
import argparse


def coordinates_list_to_BED(chrom_name: str, coordinates: list) -> str:
    '''
    function to create BED format from a list of coordinates
    takes list [[start, stop], [start, stop]]
    '''
    result = ""
    for lst in coordinates:
        result += (chrom_name + '\t' + str(lst[0]) + '\t' + str(lst[1]) + '\n')
    return result


def main():
    print('---without filter')
    coordinates = Coordinates(args.input, args.whole_genome_value,
                              args.deviation_percent)
    pseudocoordinates = coordinates.pseudocoordinates(args.coverage_column_name,
                                                      args.window_column_name,
                                                      args.repeat_window_number)
    print('---without filter')
    print(coordinates_list_to_BED(args.scaffold_name, pseudocoordinates))

    print('---filtering by distanse')
    print(coordinates_list_to_BED(args.scaffold_name, Coordinates.distanse_coordinate_filter(pseudocoordinates, args.min_region_length)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script for determining the coordinates of the pseudoautosomal region. Output of coordinates to BED file.")

    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=lambda s: metaopen(s, "rt"),
                                help="input coverage_statistics_output.csv (don`t use for STDIN)", default=stdin)

    group_additional = parser.add_argument_group('Additional options')
    group_additional.add_argument('-o', '--output', metavar='PATH', type=str, default=False,
                                  help='output file prefix')
    group_additional.add_argument('-f', '--window-size', type=int,
                                  help="the window size used in your data", default=100000)
    group_additional.add_argument('--window_column_name', type=int,
                                  help="Number of column in coverage file with window number", default=1)
    group_additional.add_argument('--coverage_column_name', type=int,
                                  help="Number of column in coverage file with mean/median coverage per window", default=2)
    group_additional.add_argument('-s', '--scaffold-name', type=str,
                                  help="Name of column in coverage file with scaffold name", default="scaffold")
    group_additional.add_argument('-m', '--whole_genome_value', type=int,
                                  help="whole genome median/mean value", default=34)
    group_additional.add_argument('-r', '--repeat_window_number', type=int,
                                  help="number of repeating windows for a given condition", default=10)
    group_additional.add_argument('-d', '--deviation_percent', type=int,
                                  help="measurement error", default=30)
    group_additional.add_argument('--min_region_length', type=int,
                                  help="minimal region length for filtration", default=15)

    args = parser.parse_args()
    main()
