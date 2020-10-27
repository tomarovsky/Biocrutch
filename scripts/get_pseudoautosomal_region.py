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


def coordinates_list_to_BED(chrom_name: str, coordinates: list) -> str:
    # takes list [[start, stop], [start, stop]]
    result = ""
    for lst in coordinates:
        result += (chrom_name + '\t' + str(lst[0]) + '\t' + str(lst[1]) + '\n')
    return result


def region_length_coordinate_filter(coordinates_list: list, min_region_length: int) -> list:
    # input list like [[start, stop], [start, stop]]
    result = []
    for lst in coordinates_list:
        if (lst[1] - lst[0]) >= min_region_length:
            result.append(lst)
    return result


def region_distanse_coordinate_filter(coordinates_list: list) -> list:
    # input list like [[start, stop], [start, stop]]
    result = []
    empty_list = True
    start_flag = True
    for lst in range(len(coordinates_list)):
        if empty_list:
            result.append(coordinates_list[lst])
            empty_list = False
            continue
        d_first = coordinates_list[lst - 1][1]
        d_second = coordinates_list[lst][0]
        distanse = d_second - d_first
        if distanse > 1:
            result.append(coordinates_list[lst])
        else:
            if start_flag == True:
                start = coordinates_list[lst - 1][0]
                start_flag = False
            stop = coordinates_list[lst][1]
            if result[-1][1] < start:
                result.append([start, stop])
                start_flag = True
            else:
                del result[-1][-1]
                result[-1].append(stop)
    return result


def main():
    deviation = args.whole_genome_value / 100 * args.deviation_percent
    minimum_coverage = args.whole_genome_value - deviation # 25.5
    maximum_coverage = args.whole_genome_value + deviation # 42.5

    coordinates = []
    repeat_window = 0
    start_coordinate = None

    for ln in args.input:
        line = ln.rstrip().split("\t")
        coverage_value = float(line[args.coverage_column_name])
        current_window = int(line[args.window_column_name])
        if coverage_value > minimum_coverage and coverage_value < maximum_coverage:
            repeat_window += 1
            if repeat_window == args.repeat_window_number and start_coordinate is None:
                start_coordinate = current_window - repeat_window + 1 #* args.window_size
                repeat_window = 0
        elif start_coordinate is not None and (coverage_value <= minimum_coverage or coverage_value >= maximum_coverage):
            stop_coordinate = current_window #* args.window_size
            coordinates.append([start_coordinate, stop_coordinate])
            start_coordinate = None
            repeat_window = 0
        else:
            repeat_window = 0


    print(coordinates)
    print('---without filter')
    print(coordinates_list_to_BED(args.scaffold_name, coordinates))
    print('---filtering by distanse')
    print(coordinates_list_to_BED(args.scaffold_name, region_distanse_coordinate_filter(coordinates)))
    print('---filtering by distanse + region length')
    print(coordinates_list_to_BED(args.scaffold_name, region_length_coordinate_filter(region_distanse_coordinate_filter(coordinates), args.min_region_length)))


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
    group_additional.add_argument('--window_column_name', type=int,
                                  help="Number of column in coverage file with window number", default=1)
    group_additional.add_argument('--coverage_column_name', type=int,
                                  help="Number of column in coverage file with mean/median coverage per window", default=2)
    group_additional.add_argument('-s', '--scaffold-name', type=str,
                                  help="Name of column in coverage file with scaffold name", default="scaffold")
    group_additional.add_argument('-m', '--whole_genome_value', type=int,
                                  help="whole genome median/mean value", default=34)
    group_additional.add_argument('-r', '--repeat_window_number', type=int,
                                  help="number of repeating windows for a given condition", default=3)
    group_additional.add_argument('-d', '--deviation_percent', type=int,
                                  help="number of repeating windows for a given condition", default=25)
    group_additional.add_argument('--min_region_length', type=int,
                                  help="minimal region length for filtration", default=15)

    args = parser.parse_args()
    main()