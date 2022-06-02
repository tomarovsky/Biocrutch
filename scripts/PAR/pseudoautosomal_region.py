#!/usr/bin/env python3
__author__ = 'tomarovsky'

from Biocrutch.Statistics.coverage_statistics.CoverageMetrics import CoveragesMetrics
from Biocrutch.Statistics.PAR.coordinator import Coordinator
from Biocrutch.Statistics.PAR.filter import Filter
from Biocrutch.Routines.routine_functions import metaopen
from collections import Counter
from sys import stdin
import argparse


def coordinates_list_to_BED(scaffold_name: str, coordinates: list) -> str:
    """
    function to create BED format from a list of coordinates
    takes list [[start, stop], [start, stop]]
    """
    result = ""
    for lst in coordinates:
        result += (scaffold_name + '\t' + str(lst[0]) + '\t' + str(lst[1]) + '\n')
    return result


def main():
    coordinates = Coordinator(args.input, float(args.whole_genome_value), args.deviation_percent, args.region_gap_size)
    coordinates_and_medians = coordinates.get_coordinates(args.window_size,
                                                  args.coverage_column_name,
                                                  args.window_column_name,
                                                  args.repeat_window_number)
    coordinates_list = coordinates_and_medians[0]
    medians_list = coordinates_and_medians[1]

    print('---- Medians between regions ---- \n', medians_list, sep="")
    print("Concatenate if the median >", round(coordinates.minimum_coverage, 2))
    print('---- Raw coordinates ---- \n', coordinates_list_to_BED(args.scaffold_name, coordinates_list), sep="")

    coordinates_merge_by_median = Filter.concat_by_median(coordinates_list, # coordinates list
                                              medians_list, # median list between regions
                                              coordinates.minimum_coverage,
                                              coordinates.maximum_coverage)
    print('---- Filtration by median ---- \n', coordinates_list_to_BED(args.scaffold_name, coordinates_merge_by_median), sep="")

    coordinates_merge_by_distance = Filter.concat_by_distance(coordinates, args.min_region_length)

    print('---- Filtration by distance ---- \n', coordinates_list_to_BED(args.scaffold_name, coordinates_merge_by_distance), sep="")

    if args.output:
        outfile = open(args.output + "_pseudoreg.bed", "w")
        outfile.writelines(coordinates_list_to_BED(args.scaffold_name, coordinates_merge_by_distance))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script for determining the coordinates of the PAR. Output of coordinates to BED file")

    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=lambda s: metaopen(s, "rt"),
                                help="input coverage_statistics_output.csv (don`t use for STDIN)", default=stdin)

    group_additional = parser.add_argument_group('Additional options')
    group_additional.add_argument('-o', '--output', metavar='PATH', type=str, default=False,
                                  help='output file prefix')
    group_additional.add_argument('-f', '--window-size', type=int,
                                  help="the window size used in your data")
    group_additional.add_argument('--window_column_name', type=int,
                                  help="number of column in coverage file with window number", default=1)
    group_additional.add_argument('--coverage_column_name', type=int,
                                  help="number of column in coverage file with mean/median coverage per window", default=2)
    group_additional.add_argument('-s', '--scaffold-name', type=str,
                                  help="name of column in coverage file with scaffold name", default="ChrX")
    group_additional.add_argument('-m', '--whole_genome_value', type=str,
                                  help="whole genome median/mean value")
    group_additional.add_argument('-r', '--repeat_window_number', type=int,
                                  help="number of repeating windows for a given condition", default=10)
    group_additional.add_argument('-g', '--region_gap_size', type=int,
                                  help="minimum allowable gap between regions for their merging (windows)", default=1)
    group_additional.add_argument('-l', '--min_region_length', type=int,
                                  help="minimum distance between regions for their merging (windows)", default=5)
    group_additional.add_argument('-d', '--deviation_percent', type=int,
                                  help="measurement error", default=30)
    args = parser.parse_args()
    main()
