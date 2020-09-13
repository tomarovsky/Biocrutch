'''
Script for calculating median, average, maximum and minimum coverage.
a. calculate stats for whole genome
b. calculate stats for each scaffold
c. calculate stats in 100 kbp and 1 Mbp stacking windows
Header-less tab-separated input file with 3 columns: scaffold_id, position(1-based), coverage.
'''
from Biocrutch.Routines.routine_functions import metaopen
from collections import Counter, OrderedDict, defaultdict
from sys import stdin
import argparse


def main():
    genome_line_counter = 0
    genome_coverage_amount = 0
    genome_coverages_amounts_dict = Counter()

    scaffold_coverages_dict = defaultdict(list)
    for line in args.input:
        # for whole genome
        genome_line_counter += 1
        line = line.rstrip().split('\t')
        genome_coverages_amounts_dict[int(line[2])] += 1
        genome_coverage_amount += int(line[2])

        # for each scaffold
        scaffold_coverages_dict[line[0]].append(int(line[2]))
        if len(scaffold_coverages_dict.keys()) > 1:
            #median
            scaffold_list_of_coverages = sorted(scaffold_coverages_dict[l])
            if len(scaffold_coverages_dict[l]) is int:
                index = len(scaffold_coverages_dict[l]) // 2
                print("scaffold", l, "median", scaffold_list_of_coverages[index])
            else:
                index_1 = int(len(scaffold_coverages_dict[l])/2 - 0.5)
                value_1 = scaffold_list_of_coverages[index_1]
                index_2 = int(len(scaffold_coverages_dict[l])/2 - 0.5)
                value_2 = scaffold_list_of_coverages[index_2]
                print("scaffold", l, "median", (value_1 + value_2) / 2)

            # average, max, min
            scaffold_line_counter = len(scaffold_coverages_dict[l])
            print("scaffold", l, "max", max(scaffold_coverages_dict[l]))
            print("scaffold", l, "min", min(scaffold_coverages_dict[l]))
            print("scaffold", l, "average", sum(scaffold_coverages_dict.pop(l))/scaffold_line_counter)
        l = line[0]

    scaffold_list_of_coverages = sorted(scaffold_coverages_dict[l])
    if len(scaffold_coverages_dict[l]) is int:
        index = int(len(scaffold_coverages_dict[l]) // 2)
        print("scaffold", l, "median", scaffold_list_of_coverages[index])
    else:
        index_1 = int(len(scaffold_coverages_dict[l])/2 - 0.5)
        value_1 = scaffold_list_of_coverages[index_1]
        index_2 = int(len(scaffold_coverages_dict[l])/2 - 0.5)
        value_2 = scaffold_list_of_coverages[index_2]
        print("scaffold", l, "median", (value_1 + value_2) / 2)

    print("scaffold", l, "max", max(scaffold_coverages_dict[l]))
    print("scaffold", l, "min", min(scaffold_coverages_dict[l]))
    scaffold_line_counter = len(scaffold_coverages_dict[l])
    print("scaffold", l, "average", sum(scaffold_coverages_dict.pop(l))/scaffold_line_counter)

    # median for whole genome
    list_of_coverages = sorted(genome_coverages_amounts_dict.keys())
    half_sum_of_coverage_amount = genome_coverage_amount/2
    count = 0
    if half_sum_of_coverage_amount % 2 != 0.0:
        for i in list_of_coverages:
            step = i * genome_coverages_amounts_dict[i]
            count += step
            if count >= half_sum_of_coverage_amount:
                genome_median = i
                break
    else:
        half_sum_of_coverage_amount += 1
        for i in list_of_coverages:
            step = [i] * genome_coverages_amounts_dict[i]
            for j in step:
                count += j
                if count == int(half_sum_of_coverage_amount):
                    genome_median = j
                    break
            step = []

    print("whole genome median", genome_median)
    print("whole genome min", list_of_coverages[0])
    print("whole genome max", list_of_coverages[-1])
    print("whole genome average", round(genome_coverage_amount/genome_line_counter, 2))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="script for calculating median, average, maximum and minimum coverage in file.bam")
    group_additional = parser.add_argument_group('Additional options')
    group_additional.add_argument('-i', '--input', type=lambda s: metaopen(s, "rt"),
                            help="input file.bam.gz (don`t use for STDIN)", default=stdin)

    args = parser.parse_args()
    main()