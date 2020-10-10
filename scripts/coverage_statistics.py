#!/usr/bin/env python3
__author__ = 'tomarovsky'
'''
Script for calculating median, average, maximum and minimum coverage.
a. calculate stats for whole genome
b. calculate stats for each scaffold
c. calculate stats in 100 kbp and 1 Mbp stacking windows
Header-less tab-separated input file with 3 columns: scaffold_id, position(1-based), coverage.
'''
from Biocrutch.Statistics.coverage_statistics.coverage_stats_functions import count_window_stats
from Biocrutch.Routines.routine_functions import metaopen
from collections import Counter
from sys import stdin
import pandas as pd
import argparse


def main():
    # created dataframe for whole genome and scaffolds stats
    df_scaffolds = pd.DataFrame(columns=['median', 'average', 'max', 'min'])
    df_frames = pd.DataFrame(columns=['scaffold', 'frame', 'median', 'average', 'max', 'min'])
    df_whole_genome = pd.DataFrame(columns=['median', 'average', 'max', 'min'])

    genome_coverages_amounts_dict = Counter()
    genome_line_counter = 0
    frame_coverages_amounts_dict = Counter()
    frame_line_counter = 0
    frame_id = -1
    index = 0
    # overlapping_frame_id = 0.5
    # overlapping_frame_coverages_amounts_dict = Counter()
    # overlapping_frame_line_counter = 0
    scaffold_coverages_dict = Counter()

    previous_scaffold_name = None

    for line in args.input:
        line = line.rstrip().split('\t')
        # for whole genome
        genome_line_counter += 1
        genome_coverages_amounts_dict[int(line[2])] += 1

        # if genome_line_counter > int(args.frame_size / 2):
        #     overlapping_frame_line_counter += 1
        #     overlapping_frame_coverages_amounts_dict[line[2]] += 1

        if previous_scaffold_name == line[0] or previous_scaffold_name == None:
            frame_line_counter += 1
            frame_coverages_amounts_dict[int(line[2])] += 1
        else:
            frame_id = -1
            frame_line_counter = 1
            frame_coverages_amounts_dict.clear()
            frame_coverages_amounts_dict[int(line[2])] += 1

        # for each scaffold
        if previous_scaffold_name != line[0] and previous_scaffold_name != None:
            df_scaffolds.loc[previous_scaffold_name] = count_window_stats(scaffold_coverages_dict)
            scaffold_coverages_dict.clear()

        # for window (non-overlapping)
        if frame_line_counter == args.frame_size:
            index += 1
            frame_id += 1
            metrics = [previous_scaffold_name, frame_id] + count_window_stats(frame_coverages_amounts_dict)
            df_frames.loc[index] = metrics
            frame_coverages_amounts_dict.clear()
            frame_line_counter = 0

        # # for window (overlapping)
        # elif overlapping_frame_line_counter == args.frame_size:
        #     overlapping_frame_id += 1
        #     metrics = [previous_scaffold_name] + count_window_stats(overlapping_frame_coverages_amounts_dict)
        #     df_frames.loc[overlapping_frame_id] = metrics
        #     overlapping_frame_coverages_amounts_dict.clear()
        #     overlapping_frame_line_counter = 0

        scaffold_coverages_dict[int(line[2])] += 1
        previous_scaffold_name = line[0]

    # processing residual data after a cycle for each scaffold and whole genome
    df_scaffolds.loc[previous_scaffold_name] = count_window_stats(scaffold_coverages_dict)
    df_whole_genome.loc['whole_genome'] = count_window_stats(genome_coverages_amounts_dict)

    #for print dataframe to terminal
    print (df_scaffolds.rename_axis('scaffold').reset_index())
    print (df_frames)
    print (df_whole_genome.rename_axis('genome').reset_index())

    if args.output: # create a report.csv
        df_scaffolds.rename_axis('scaffold').reset_index().to_csv(args.output + "_scaffolds_stats.csv", encoding='utf-8', sep='\t')
        df_frames.to_csv(args.output + '_' + str(args.frame_size) + "_windows_stats.csv", encoding='utf-8', sep='\t')
        df_whole_genome.rename_axis('genome').reset_index().to_csv(args.output + "_whole_genome_stats.csv", encoding='utf-8', sep='\t')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="script for calculating median, average, maximum and minimum coverage in file.tab. Report to files.csv")

    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=lambda s: metaopen(s, "rt"),
                                  help="input file.bam.gz (don`t use for STDIN)", default=stdin)

    group_additional = parser.add_argument_group('Additional options')
    group_additional.add_argument('-f', '--frame-size', type=int,
                                  help="calculate stats in 100 kbp and 1 Mbp stacking windows", default=500)
    group_additional.add_argument('-o', '--output', metavar='PATH', type=str,
                                  help='output file prefix without frame size')
    args = parser.parse_args()
    main()
