#!/usr/bin/env python3
from collections import defaultdict
from Bio import AlignIO
from sys import stdin
import argparse


def main():
    with open(args.output, "a") as outfile:
        AlignIO.convert(args.input, "fasta", outfile, args.format, args.type)
        outfile.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="script to convert FASTA into some format")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str, default=stdin, help="input concat FASTA file or stdin")
    group_required.add_argument('-f', '--format', type=str,
                                help="output file format ('nexus', 'phylip' or any other supported by Bio.AlignIO)")
    group_required.add_argument('-t', '--type', type=str, help="molecular type (DNA, RNA or protein)")
    group_required.add_argument('-o', '--output', type=str, help="output file name")
    args = parser.parse_args()
    main()