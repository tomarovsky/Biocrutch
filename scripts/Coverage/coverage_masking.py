#!/usr/bin/env python3
__author__ = 'tomarovsky'
'''
Get coverage mask. You can use pipeline:
zcat SAMPLE.per-base.bed.gz | awk '{if( ($4 > 2.5 * whole_genome_median) || ($4 < 0.5 * whole_genome_median) print $0)}' | covarege_masking.py
WHOLE_MEDIAN=$(cat *_whole_genome_stats.csv | sed -n 2p | awk '{print $2}');
zcat SAMPLE.per-base.bed.gz | awk '{if( ($4 > 2.5 * $WHOLE_MEDIAN) || $4 < 0.5 * $WHOLE_MEDIAN) print $0)}' | covarege_masking.py
'''
from Biocrutch.Routines.routine_functions import metaopen, metaoutput
import argparse
from sys import stdin


def main():
    outfile = metaopen(metaoutput(args.output, ".bed.gz"), "wt")
    scaffold_name = None
    previous_scaffold_name = None
    global_start = None
    start = None
    stop = None

    for i in args.input:
        line = i.strip().split()
        scaffold_name = line[0]
        if scaffold_name == previous_scaffold_name or scaffold_name is None:
            if int(line[3]) > (2.5 * args.whole_median) or int(line[3]) < (0.5 * args.whole_median):
                if global_start is None:
                    global_start = line[1]
                    stop = line[2]
                    continue
                else:
                    start = line[1]
                    # print (stop, start)
                    if stop == start:
                        stop = line[2]
                    else:
                        outline = "\t".join([scaffold_name, global_start, stop])
                        outfile.write(outline)
                        global_start = line[1]
                        stop = line[2]
                        # print (outline)
        else:
            if previous_scaffold_name is not None:
                outline = "\t".join([previous_scaffold_name, global_start, stop])
                # print (outline)
            global_start = line[1]
            stop = line[2]
            start = None
        previous_scaffold_name = line[0]
    outline = "\t".join([previous_scaffold_name, global_start, stop])
    # print (outline)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="get the coverage mask through the pipeline")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=lambda s: metaopen(s, "rt"),
                                help="input coverage_statistics_output.csv (don`t use for STDIN)", default=stdin)
    group_required.add_argument('-w', '--whole-median', type=int,
                                help="whole median value")
    group_required.add_argument('-o', '--output', type=str,
                                help="output file prefix")
    args = parser.parse_args()
    main()