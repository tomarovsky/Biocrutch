#!/usr/bin/env python3
__author__ = 'tomarovsky'
'''
Script for calculating median, average, maximum and minimum coverage.
a. calculate stats for whole genome
b. calculate stats for each scaffold
c. calculate stats in 100 kbp and 1 Mbp stacking windows
Header-less tab-separated input file with 3 columns: scaffold_id, position(1-based), coverage.
'''
from Biocrutch.Routines.routine_functions import metaopen
from collections import Counter
from sys import stdin
import pandas as pd
import argparse


def window_stats(coverages_amounts_dict: Counter) -> list:
    sum_of_coverages = 0
    for key, value in coverages_amounts_dict.items():
        sum_of_coverages += key * value

    keys_coverages = sorted(coverages_amounts_dict.keys())
    sum_values_coverages = sum(coverages_amounts_dict.values())
    half_sum_values_coverages = sum_values_coverages / 2
    count = 0

    if sum_values_coverages % 2 != 0:
        for i in range(len(keys_coverages)):
            count += coverages_amounts_dict[keys_coverages[i]]
            if count >= half_sum_values_coverages:
                genome_median = keys_coverages[i]
                return [genome_median, round(sum_of_coverages/sum(coverages_amounts_dict.values()), 2),
                        keys_coverages[-1],
                        keys_coverages[0]]
    else:
        half_sum_values_coverages = int(half_sum_values_coverages)
        for i in range(len(keys_coverages)):
            count += coverages_amounts_dict[keys_coverages[i]]
            if count == half_sum_values_coverages:
                genome_median = (keys_coverages[i] + keys_coverages[i+1])/2
                return [genome_median, round(sum_of_coverages/sum(coverages_amounts_dict.values()), 2),
                        keys_coverages[-1], keys_coverages[0]]
            elif count > half_sum_values_coverages:
                genome_median = keys_coverages[i]
                return [genome_median, round(sum_of_coverages/sum(coverages_amounts_dict.values()), 2),
                        keys_coverages[-1],
                        keys_coverages[0]]


def pretty_printer(dataframe):
    # print and retyping for values
    dataframe[['max', 'min']] = dataframe[['max', 'min']].astype(int)
    dataframe['median'] = dataframe['median'].astype(str)
    dataframe.average = dataframe.average.round(2)
    print (dataframe)
    return dataframe


def main():
    # created dataframe for whole genome and scaffolds stats
    df_scaffolds = pd.DataFrame(columns=['median', 'average', 'max', 'min'])
    df_frames = pd.DataFrame(columns=['median', 'average', 'max', 'min'])
    df_whole_genome = pd.DataFrame(columns=['median', 'average', 'max', 'min'])

    genome_coverages_amounts_dict = Counter()
    genome_line_counter = 0
    frame_coverages_amounts_dict = Counter()
    frame_line_counter = 0
    frame_id = 0
    overlapping_frame_id = 0.5
    overlapping_frame_coverages_amounts_dict = Counter()
    overlapping_frame_line_counter = 0
    scaffold_coverages_dict = Counter()

    for line in args.input:
        # for whole genome
        line = line.rstrip().split('\t')
        genome_line_counter += 1
        frame_line_counter += 1
        frame_coverages_amounts_dict[int(line[2])] += 1
        scaffold_coverages_dict[int(line[2])] += 1
        genome_coverages_amounts_dict[int(line[2])] += 1
        if genome_line_counter >= int(args.frame_size / 2): # !!!!!!!add args
            overlapping_frame_line_counter += 1
            overlapping_frame_coverages_amounts_dict[int(line[2])] += 1

        # for each scaffold
        try:
            if previous_scaffold_name != line[0]:
                df_scaffolds.loc[previous_scaffold_name] = window_stats(scaffold_coverages_dict)
                scaffold_coverages_dict.clear()
        except UnboundLocalError:
            pass

        # for window (non-overlapping)
        if frame_line_counter == args.frame_size:
            frame_id += 1
            df_frames.loc['frame_'+str(frame_id)] = window_stats(frame_coverages_amounts_dict)
            frame_coverages_amounts_dict.clear()
            frame_line_counter = 0

        # for window (overlapping)
        elif overlapping_frame_line_counter == args.frame_size:
            overlapping_frame_id += 1
            df_frames.loc['frame_'+str(overlapping_frame_id)] = window_stats(overlapping_frame_coverages_amounts_dict)
            overlapping_frame_coverages_amounts_dict.clear()
            overlapping_frame_line_counter = 0

        previous_scaffold_name = line[0]

    # processing residual data after a cycle for each scaffold and whole genome
    df_scaffolds.loc[previous_scaffold_name] = window_stats(scaffold_coverages_dict)
    df_whole_genome.loc['whole_genome'] = window_stats(genome_coverages_amounts_dict)

    #for print dataframe to terminal
    df_scaffolds = pretty_printer(df_scaffolds)
    df_frames = pretty_printer(df_frames)
    df_whole_genome = pretty_printer(df_whole_genome)

    if args.output: # create a report.csv
        df_scaffolds.rename_axis('scaffold').reset_index().to_csv(args.output + "_scaffolds_stats.csv", encoding='utf-8', sep='\t')
        df_frames.rename_axis('frame').reset_index().to_csv(args.output + "_windows_stats.csv", encoding='utf-8', sep='\t')
        df_whole_genome.rename_axis('genome').reset_index().to_csv(args.output + "_whole_genome_stats.csv", encoding='utf-8', sep='\t')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="script for calculating median, average, maximum and minimum coverage in file.tab. Report to files.csv")

    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=lambda s: metaopen(s, "rt"),
                                  help="input file.bam.gz (don`t use for STDIN)", default=stdin)

    group_additional = parser.add_argument_group('Additional options')
    group_additional.add_argument('-f', '--frame-size', type=int,
                                  help="calculate stats in 100 kbp and 1 Mbp stacking windows", default=500)
    group_additional.add_argument('-o', '--output', metavar='PATH', type=str,
                                  help='output file prefix')
    args = parser.parse_args()
    main()
