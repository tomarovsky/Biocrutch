#!/usr/bin/env python3
__author__ = 'tomarovsky'
from argparse import ArgumentParser


CSCAFFOLDS = ['HiC_scaffold_1', 'HiC_scaffold_2', 'HiC_scaffold_3',
'HiC_scaffold_4', 'HiC_scaffold_5', 'HiC_scaffold_6', 'HiC_scaffold_7',
'HiC_scaffold_8', 'HiC_scaffold_9', 'HiC_scaffold_10', 'HiC_scaffold_11',
'HiC_scaffold_12', 'HiC_scaffold_13', 'HiC_scaffold_14', 'HiC_scaffold_15',
'HiC_scaffold_16', 'HiC_scaffold_17', 'HiC_scaffold_18', 'HiC_scaffold_19']


def lengths_to_dict(lenfile):
    lengths_dict = {}
    with open(lenfile, 'r') as lengths:
        for line in lengths:
            line = line.strip().split("\t")
            lengths_dict[line[0]] = int(line[1])
    return lengths_dict


def main():
    lengths_dict = lengths_to_dict(args.len_file)
    result = []
    with open(args.dups_file, 'r') as dups:
        for line in dups:
            line = line.strip().split("\t")
            if len(line) == 5 and line[3] == "HAPLOTIG":
                if line[0] in CSCAFFOLDS and line[4] not in CSCAFFOLDS:
                    result.append(line[4])
                if line[4] in CSCAFFOLDS and line[0] not in CSCAFFOLDS:
                    result.append(line[0])
                if line[0] not in CSCAFFOLDS and line[4] not in CSCAFFOLDS:
                    if lengths_dict[line[0]] == lengths_dict[line[4]]:
                        result.append(line[4])
                    elif lengths_dict[line[0]] > lengths_dict[line[4]]:
                        result.append(line[0])
                    else:
                        result.append(line[4])
    result = list(set(result))
    print("removed.ids file contains: ", len(result))
    with open(args.output, 'w') as outfile:
        outfile.write("\n".join(result))


if __name__ == "__main__":
    parser = ArgumentParser(description="FASTA filtering from haplotigs based on dups.bed")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-l', '--len-file', type=str, help="FASTA.len file")
    group_required.add_argument('-d', '--dups-file', type=str, help="dups.bed file (from purge_dups output)")
    group_required.add_argument('-o', '--output', type=str, help="removed.ids file")
    args = parser.parse_args()
    main()