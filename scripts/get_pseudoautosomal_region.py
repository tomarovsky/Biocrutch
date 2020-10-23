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

def coordinate_filter(coordinates_list: list) -> list:
    # function for filtering coordinates.
    # Accepts most likely a list [[start, stop], [start, stop]]
    start = None
    stop = None
    result = []
    for lst in range(len(coordinates_list)):
        if lst >= 1:
            stop = coordinates_list[lst - 1][1]
            start = coordinates_list[lst][0]
            distanse = start - stop
            if distanse > 10:
                if coordinates_list[lst - 1] not in result:
                    result.append(coordinates_list[lst - 1])
                result.append(coordinates_list[lst])
    return result


def coordinates_list_to_BED(chrom_name: str, coordinates_list: list) -> str:
    # function to convert coordinates to BED format.
    # Accepts most likely a list [[start, stop], [start, stop]]
    result = ""
    for lst in coordinates_list:
        result += (chrom_name + '\t' + str(lst[0]) + '\t' + str(lst[1]) + '\n')
    return result


def main():
    deviation = args.whole_genome_value / 100 * args.deviation_percent
    minimum_coverage = args.whole_genome_value - deviation # 25.5
    maximum_coverage = args.whole_genome_value + deviation # 42.5

    coordinates = []
    repeat_window = 0
    start_coordinate = None

    for line in args.input:
        line = line.rstrip().split("\t")
        coverage_value = float(line[args.coverage_column_name])
        current_window = int(line[args.window_column_name])
        if coverage_value > minimum_coverage and coverage_value < maximum_coverage:
            repeat_window += 1
            if repeat_window == args.repeat_window_number and start_coordinate is None:
                start = current_window - repeat_window + 1
                start_coordinate = start #* args.window_size
                repeat_window = 0
        elif start_coordinate is not None:
            stop = current_window
            stop_coordinate = stop #* args.window_size
            coordinates.append([start_coordinate, stop_coordinate])
            start_coordinate = None
        else:
            repeat_window = 0


    print(coordinates)
    coord = coordinate_filter(coordinates)
    print(coordinates_list_to_BED(args.scaffold_name, coord))
    print('------------')
    print(coordinates_list_to_BED(args.scaffold_name, coordinates))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script for determining the coordinates of the pseudoautosomal region. Output of coordinates to BED file.")

    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=lambda s: metaopen(s, "rt"),
                                  help="input file.bam.gz (don`t use for STDIN)", default=stdin)

    group_additional = parser.add_argument_group('Additional options')
    group_additional.add_argument('-o', '--output', metavar='PATH', type=str, default=False,
                                  help='output file prefix')
    group_additional.add_argument('-f', '--window-size', type=int,
                                  help="the window size used in your data", default=100000)
    group_additional.add_argument('--coverage_column_name', type=int,
                                  help="Number of column in coverage file with mean/median coverage per window", default=2)
    group_additional.add_argument('--window_column_name', type=int,
                                  help="Number of column in coverage file with window number", default=1)
    group_additional.add_argument('-s', '--scaffold-name', type=str,
                                  help="Name of column in coverage file with scaffold name", default="scaffold")
    group_additional.add_argument('-m', '--whole_genome_value', type=int,
                                  help="whole genome median/mean value", default=34)
    group_additional.add_argument('-r', '--repeat_window_number', type=int,
                                  help="number of repeating windows for a given condition", default=3)
    group_additional.add_argument('-d', '--deviation_percent', type=int,
                                  help="number of repeating windows for a given condition", default=25)

    args = parser.parse_args()
    main()