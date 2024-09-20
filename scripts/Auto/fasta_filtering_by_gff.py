#!/usr/bin/env python3
__author__ = "tomarovsky"
from argparse import ArgumentParser

from Bio import SeqIO


def extract_scaffolds_from_gff(gff_file):
    scaffolds = set()
    with open(gff_file, "r") as gff:
        for line in gff:
            if line.startswith("#"):
                continue
            columns = line.split("\t")
            scaffold = columns[0]
            scaffolds.add(scaffold)
    return scaffolds


def filter_fasta_by_scaffolds(fasta_file, scaffolds, output_fasta):
    with open(output_fasta, "w") as output_handle:
        fasta_sequences = SeqIO.parse(fasta_file, "fasta")
        filtered_sequences = (record for record in fasta_sequences if record.id in scaffolds)
        SeqIO.write(filtered_sequences, output_handle, "fasta")


def main():
    scaffolds = extract_scaffolds_from_gff(args.gff)
    filter_fasta_by_scaffolds(args.fasta, scaffolds, args.output)


if __name__ == "__main__":
    parser = ArgumentParser(description="FASTA filtering by ids from tabular file")
    group_required = parser.add_argument_group("Required options")
    group_required.add_argument("-f", "--fasta", type=str, help="FASTA file")
    group_required.add_argument("-g", "--gff", type=str, help="GFF file")
    group_required.add_argument("-o", "--output", type=str, help="output FASTA name")
    args = parser.parse_args()
    main()
