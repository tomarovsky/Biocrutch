#!/usr/bin/env python3
__author__ = 'tomarovsky'
from argparse import ArgumentParser

def reverselist_to_list(reverse_file):
    reverselist = []
    with open(reverse_file, 'r') as revs:
        for line in revs:
            reverselist.append(line.strip())
    return reverselist

def lengths_to_dict(lenfile):
    lengths_dict = {}
    with open(lenfile, 'r') as lengths:
        for line in lengths:
            line = line.strip().split("\t")
            lengths_dict[line[0]] = int(line[1])
    return lengths_dict

def main():
    reverselist = reverselist_to_list(args.reverse_file)
    lengths_dict = lengths_to_dict(args.len_file)
    outfile = open(args.output, 'wt')
    count = 0
    for scaffold, length in lengths_dict.items():
        count += 1
        if scaffold in reverselist:
            chain = f"chain 1000000 {scaffold}_rc {length} - 0 {length} {scaffold} {length} + 0 {length} {count}\n{length}\n"
        else:
            chain = f"chain 1000000 {scaffold} {length} + 0 {length} {scaffold} {length} + 0 {length} {count}\n{length}\n"
        outfile.write(chain)


if __name__ == "__main__":
    parser = ArgumentParser(description="Generate a Chain file for reverse FASTA")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--len-file', type=str, help="FASTA.len file")
    group_required.add_argument('-r', '--reverse-file', type=str, help=".reverselist file")
    group_required.add_argument('-o', '--output', type=str, help="Chain file with scaffolds from *.reverselist file")
    args = parser.parse_args()
    main()