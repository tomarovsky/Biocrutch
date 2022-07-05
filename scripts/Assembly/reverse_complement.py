#!/usr/bin/env python3
__author__ = 'tomarovsky'
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Seq import Seq


def reverselist_to_list(reverse_file):
    reverselist = []
    with open(reverse_file, 'r') as revs:
        for line in revs:
            reverselist.append(line.strip())
    return reverselist

def main():
    reverselist = reverselist_to_list(args.reverse_file)
    fasta_sequences = SeqIO.parse(open(args.input),'fasta')
    reverse_fasta = open(args.output, "wt")
    for scaffold in fasta_sequences:
        if scaffold.id in reverselist:
            scaffold.seq = scaffold.seq.reverse_complement()
            scaffold.id += "_rc"
            scaffold.description += "_rc"
        SeqIO.write(scaffold, reverse_fasta, "fasta")


if __name__ == "__main__":
    parser = ArgumentParser(description="Reverse complement same scaffolds in FASTA")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str, help="FASTA file")
    group_required.add_argument('-r', '--reverse-file', type=str, help=".reverselist file")
    group_required.add_argument('-o', '--output', type=str, help="FASTA with reverse complement scaffolds")
    args = parser.parse_args()
    main()