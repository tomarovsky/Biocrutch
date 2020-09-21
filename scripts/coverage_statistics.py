'''
Script for calculating median, average, maximum and minimum coverage.
a. calculate stats for whole genome
b. calculate stats for each scaffold
c. calculate stats in 100 kbp and 1 Mbp stacking windows
Header-less tab-separated input file with 3 columns: scaffold_id, position(1-based), coverage.
'''
from Biocrutch.Routines.routine_functions import metaopen
from collections import Counter, OrderedDict, defaultdict
import pandas as pd
from sys import stdin
import argparse


def each_scaffold_stats(scaffold_coverages_dict, l):
    # use to calculate stats for each scaffold
    tmp_lst_to_df = []
    # median
    scaffold_list_of_coverages = sorted(scaffold_coverages_dict[l])
    if len(scaffold_coverages_dict[l]) is int:
        index = len(scaffold_coverages_dict[l]) // 2
        tmp_lst_to_df.append(scaffold_list_of_coverages[index])
    else:
        index_1 = int(len(scaffold_coverages_dict[l])/2 - 0.5)
        value_1 = scaffold_list_of_coverages[index_1]
        index_2 = int(len(scaffold_coverages_dict[l])/2 - 0.5)
        value_2 = scaffold_list_of_coverages[index_2]
        tmp_lst_to_df.append((value_1 + value_2) / 2)
    # average, max, min
    scaffold_line_counter = len(scaffold_coverages_dict[l])
    tmp_lst_to_df.append(sum(scaffold_coverages_dict[l])/scaffold_line_counter)
    tmp_lst_to_df.append(max(scaffold_coverages_dict[l]))
    # .pop() to clear the dictionary
    tmp_lst_to_df.append(min(scaffold_coverages_dict.pop(l)))
    return tmp_lst_to_df


def frame_stats(coverages_amounts_dict, coverage_amount, line_counter):
    # use to calculate stats for whole genome and stacking windows
    # median
    list_of_coverages = sorted(coverages_amounts_dict.keys())
    half_sum_of_coverage_amount = coverage_amount/2
    count = 0

    for i in range(len(list_of_coverages)):
        for j in range(coverages_amounts_dict[list_of_coverages[i]]):
            count += list_of_coverages[i]
            if count > half_sum_of_coverage_amount:
                if j != 0:
                    genome_median = list_of_coverages[i]
                    return [genome_median,
                            round(coverage_amount/line_counter, 2),
                            list_of_coverages[-1],
                            list_of_coverages[0]]
                else:
                    genome_median = (list_of_coverages[i] + list_of_coverages[i - 1]) / 2
                    return [genome_median,
                            round(coverage_amount/line_counter, 2),
                            list_of_coverages[-1],
                            list_of_coverages[0]]


def main():
    # created dataframe for whole genome and scaffolds stats
    df_whole_and_scaffolds = pd.DataFrame(
        columns=['median', 'average', 'max', 'min'])
    df_stacking_windows = pd.DataFrame(
        columns=['median', 'average', 'max', 'min'])

    genome_coverages_amounts_dict = Counter()
    genome_coverage_amount = 0
    genome_line_counter = 0
    frame_coverages_amounts_dict = Counter()
    frame_coverage_amount = 0
    frame_line_counter = 0
    frame_id = 0
    overlapping_frame_id = 0.5
    overlapping_frame_coverages_amounts_dict = Counter()
    overlapping_frame_line_counter = 0
    overlapping_frame_coverage_amount = 0
    scaffold_coverages_dict = defaultdict(list)

    for line in args.input:
        # for whole genome
        genome_line_counter += 1
        frame_line_counter += 1
        line = line.rstrip().split('\t')
        frame_coverages_amounts_dict[int(line[2])] += 1
        frame_coverage_amount += int(line[2])
        genome_coverages_amounts_dict[int(line[2])] += 1
        genome_coverage_amount += int(line[2])
        if genome_line_counter >= int(args.frame_size / 2):
            overlapping_frame_line_counter += 1
            overlapping_frame_coverage_amount += int(line[2])
            overlapping_frame_coverages_amounts_dict[int(line[2])] += 1

        # for each scaffold
        scaffold_coverages_dict[line[0]].append(int(line[2]))
        if len(scaffold_coverages_dict.keys()) > 1:  # for each scaffold
            df_whole_and_scaffolds.loc[l] = each_scaffold_stats(
                scaffold_coverages_dict, l)

        # for each stacking window (non-overlapping)
        if frame_line_counter == args.frame_size:
            frame_id += 1
            df_stacking_windows.loc['frame_'+str(frame_id)] = frame_stats(
                frame_coverages_amounts_dict, frame_coverage_amount, frame_line_counter)
            frame_coverages_amounts_dict = Counter()
            frame_coverage_amount = 0
            frame_line_counter = 0

        # for each stacking window (overlapping)
        elif overlapping_frame_line_counter == args.frame_size:
            overlapping_frame_id += 1
            df_stacking_windows.loc['frame_'+str(overlapping_frame_id)] = frame_stats(
                overlapping_frame_coverages_amounts_dict, overlapping_frame_coverage_amount, overlapping_frame_line_counter)
            overlapping_frame_coverages_amounts_dict = Counter()
            overlapping_frame_coverage_amount = 0
            overlapping_frame_line_counter = 0

        l = line[0]

    # for each scaffold
    df_whole_and_scaffolds.loc[l] = each_scaffold_stats(
        scaffold_coverages_dict, l)
    # whole genome stats to df
    df_whole_and_scaffolds.loc['whole_genome'] = frame_stats(
        genome_coverages_amounts_dict, genome_coverage_amount, genome_line_counter)

    df_whole_and_scaffolds[['max', 'min']] = df_whole_and_scaffolds[['max', 'min']].astype(int)
    df_whole_and_scaffolds['median'] = df_whole_and_scaffolds['median'].astype(str)
    df_whole_and_scaffolds.average = df_whole_and_scaffolds.average.round(2)
    df_stacking_windows[['max', 'min']] = df_stacking_windows[['max', 'min']].astype(int)
    df_stacking_windows['median'] = df_stacking_windows['median'].astype(str)
    df_stacking_windows.average = df_stacking_windows.average.round(2)

    print(df_whole_and_scaffolds)
    print(df_stacking_windows)

    if args.output:  # create a report.csv
        df_whole_and_scaffolds.to_csv(
            args.output + "_whole_and_scaffolds.csv", encoding='utf-8')
        df_stacking_windows.to_csv(
            args.output + "_stacking_windows.csv", encoding='utf-8')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="script for calculating median, average, maximum and minimum coverage in file.bam. Report to files.csv")
    group_additional = parser.add_argument_group('Additional options')
    group_additional.add_argument('-i', '--input', type=lambda s: metaopen(s, "rt"),
                                  help="input file.bam.gz (don`t use for STDIN)", default=stdin)
    group_additional.add_argument('-f', '--frame-size', type=int,
                                  help="calculate stats in 100 kbp and 1 Mbp stacking windows", default=500)
    group_additional.add_argument('-o', '--output', metavar='PATH', type=str,
                                  help='output file prefix')
    args = parser.parse_args()
    main()
