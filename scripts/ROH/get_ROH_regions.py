#!/usr/bin/env python3
__author__ = 'tomarovsky'
from argparse import ArgumentParser
import pandas as pd


def filter_by_level_of_heterozygosity(df, threshold, max_threshold_factor, max_adjacent_windows):
    # step 1: set max threshold
    df = df.loc[df.iloc[:, -1] <= (args.threshold * args.max_threshold_factor)].reset_index(drop=True)
    # step 2: remove adjacent windows exceeding the threshold if there are more than max_adjacent_windows
    remove_list = []
    count = 0
    for i in range(len(df)):
        if df.iloc[i][-1] >= args.threshold:
            count += 1
        else:
            if count >= args.max_adjacent_windows:
                for j in range(i - count, i):
                    remove_list.append(j)
            count = 0
    df.drop(df.index[remove_list], inplace=True)
    df = df.reset_index(drop=True)
    return df


def concat_regions_by_distance(df, min_distance):
    '''
    merge coordinates that are less than the specified distance + calculation of ROH lengths
    '''
    prev_scaffold = None
    result = []
    start_flag = True
    for scaffold in df['scaffold'].unique():
        for i in range(len(df)):
            if df.iloc[i]["scaffold"] == scaffold:
                if prev_scaffold != scaffold:
                    result.append([scaffold, df.iloc[i]["start"], df.iloc[i]["end"]])
                    prev_scaffold = scaffold
                    continue
                distance = df.iloc[i]["start"] - df.iloc[i-1]["end"]
                if distance > min_distance:
                    result.append([scaffold, df.iloc[i]["start"], df.iloc[i]["end"]])
                else:
                    del result[-1][-1]
                    result[-1].append(df.iloc[i]["end"])
    roh_df = pd.DataFrame(result, columns=['scaffold', 'start', 'end'])
    roh_df['length'] = roh_df['end'] - roh_df['start']
    print(roh_df)
    return roh_df


def main():
    df = pd.read_csv(args.input, sep="\t")
    df = filter_by_level_of_heterozygosity(df, args.threshold, args.max_threshold_factor, args.max_adjacent_windows)
    print(df[df['scaffold'] == 'chr12'].to_string())
    roh_df = concat_regions_by_distance(df, args.min_distance)
    roh_df.to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    parser = ArgumentParser(description="script for obtaining ROH regions from variant window densities from MACE BED files")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str, help="MACE BED file")
    group_required.add_argument('-o', '--output', type=str, help="outfile name")
    group_additional = parser.add_argument_group('Additional options')
    group_additional.add_argument('-t', '--threshold', type=int, default=0.05)
    group_additional.add_argument('-m', '--max_threshold_factor', type=int, default=2)
    group_additional.add_argument('-a', '--max_adjacent_windows', type=int, default=5)
    group_additional.add_argument('-d', '--min_distance', type=int, default=500000)
    args = parser.parse_args()
    main()
