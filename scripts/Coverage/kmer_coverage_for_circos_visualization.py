#!/usr/bin/env python3
__author__ = 'tomarovsky'

from Biocrutch.Routines.routine_functions import metaopen, metaoutput
import argparse
from sys import stdin


def main():
    outfile = metaopen(metaoutput(args.output, ".csv"), "wt")
    for line in args.input:
        scaffold_name, start, stop, coverage_value = args.input.readline().strip().split()
        coverage_value = int(coverage_value)
        if (coverage_value > (2.5 * args.whole_median)) or coverage_value < (0.5 * args.whole_median):
            break

    for i in args.input:
        line = i.strip().split()
        line[3] = int(line[3])
        if (line[3] > (2.5 * args.whole_median)) or (line[3] < (0.5 * args.whole_median)):
            if line[0] == scaffold_name:
                if stop == line[1]:
                    stop = line[2]
                    continue
                else:
                    outfile.write("\t".join([scaffold_name, start, stop]) + "\n")
                    start = line[1]
                    stop = line[2]
            else:
                outfile.write("\t".join([scaffold_name, start, stop]) + "\n")
                scaffold_name, start, stop = line[:3]
    outfile.write("\t".join([scaffold_name, start, stop]) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="get kmers coverage statistics for circos")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=lambda s: metaopen(s, "rt"),
                                help="input splited_fasta_1.fasta", default=stdin)
    group_required.add_argument('-o', '--output', type=str,
                                help="output file prefix")
    args = parser.parse_args()
    main()