#!/usr/bin/env python3
__author__ = 'tomarovsky'

import argparse
from sys import stdin
from collections import Counter
from Biocrutch.Routines.routine_functions import metaopen, metaoutput
from Biocrutch.Statistics.coverage_statistics.CoverageMetrics import CoveragesMetrics


def main():
    outfile = metaopen(metaoutput(args.output, ".csv"), "wt")
    frame_coverages_amounts_dict = Counter()
    line_counter = 0
    frame_line_counter = 0

    for ln in args.input:
        line = ln.strip().split()
        line_counter += 1
        frame_line_counter += 1
        frame_coverages_amounts_dict[int(line[1])] += 1
        if frame_line_counter == args.frame_size:
            start = line_counter - args.frame_size + 1
            stop = line_counter
            metrics = CoveragesMetrics(frame_coverages_amounts_dict)
            coverage = metrics.median_value()
            outfile.write("\t".join(["MT", str(start), str(stop), str(coverage)]) + "\n")
            frame_coverages_amounts_dict.clear()
            frame_line_counter = 0
    if frame_coverages_amounts_dict:
        start = line_counter - sum(frame_coverages_amounts_dict.values()) + 1
        stop = line_counter
        metrics = CoveragesMetrics(frame_coverages_amounts_dict)
        coverage = metrics.median_value()
        outfile.write("\t".join(["MT", str(start), str(stop), str(coverage)]) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="get kmers coverage statistics for circos")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=lambda s: metaopen(s, "rt"),
                                help="input splited_fasta_1.fasta", default=stdin)
    group_required.add_argument('-k', '--kmer-length', type=int,
                                help="kmer length", default=23)
    group_required.add_argument('-f', '--frame-size', type=int,
                                help="non-overlapping frame size", default=10)
    group_required.add_argument('-o', '--output', type=str,
                                help="output file prefix")
    args = parser.parse_args()
    main()