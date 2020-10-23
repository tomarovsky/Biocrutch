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

def coordinates_list_to_BED(chrom_name: str, coordinates_list: list) -> str:
    # function to convert output to BED format. Accepts most likely a list [[start, stop], [start, stop]]
    result = ""
    for lst in coordinates_list:
        result += (chrom_name + '\t' + str(lst[0]) + '\t' + str(lst[1]) + '\n')
    return result


def main():
    deviation_percent = 25 # the percent of deviation
    deviation = args.whole_genome_value / 100 * deviation_percent
    minimum_coverage = args.whole_genome_value - deviation # 25.5
    maximum_coverage = args.whole_genome_value + deviation # 42.5

    coordinates = []
    repeat_frame = 0
    start_coordinate = None

    for line in args.input:
        line = line.rstrip().split("\t")
        coverage_value = float(line[args.coverage_column_name])
        current_frame = int(line[args.window_column_name])
        if coverage_value > minimum_coverage and coverage_value < maximum_coverage:
            repeat_frame += 1
            if repeat_frame == args.repeat_frame_number and start_coordinate is None:
                start = current_frame - repeat_frame + 1
                start_coordinate = start #* args.frame_size
                repeat_frame = 0
        elif start_coordinate is not None:
            stop = current_frame
            stop_coordinate = stop #* args.frame_size
            coordinates.append([start_coordinate, stop_coordinate])
            start_coordinate = None
        else:
            repeat_frame = 0


    print(coordinates)
    print(coordinates_list_to_BED(args.scaffold_name, coordinates))
        







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
    group_additional.add_argument('--window_column_name', type=int,
                                  help="Number of column in coverage file with mean/median coverage per window", default=2)
    group_additional.add_argument('-s', '--scaffold-name', type=str,
                                  help="Name of column in coverage file with scaffold name", default="scaffold")
    group_additional.add_argument('-m', '--whole_genome_value', type=int,
                                  help="whole genome median/mean value", default=34)
    group_additional.add_argument('-r', '--repeat_frame_number', type=int,
                                  help="number of repeating windows for a given condition", default=3)

    args = parser.parse_args()
    main()