#!/usr/bin/env python3
__author__ = 'tomarovsky'
from argparse import ArgumentParser

def main():
    pass

if __name__ == "__main__":
    parser = ArgumentParser(description="Reverse complement same scaffolds in FASTA")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str, help="FASTA file")
    group_required.add_argument('-r', '--reverse-file', type=str, help=".reverselist file")
    group_required.add_argument('-o', '--output', type=str, help="FASTA with reverse complement scaffolds")
    args = parser.parse_args()
    main()