#!/usr/bin/env python3
__author__ = 'tomarovsky'
from argparse import ArgumentParser
from Biocrutch.Parsers.fasta_opener import Fasta_opener
from Biocrutch.Routines.routine_functions import metaopen, metaoutput


def main():
    outfile = metaopen(metaoutput(args.output, '.faa'), 'wt')
    csv = metaopen(args.tabular_file, 'rt').readlines()[1:]
    data = Fasta_opener(args.input)
    faa_dict = data.parse_sequences_without_join()

    geneid_prev = None
    count = 0
    for ln in csv:
        line = ln.strip().split(',')
        geneid, length = line[5], int(line[9])
        if geneid_prev is None:
            geneid_prev = geneid
            count, protein = length, line[8][1:-1]
            continue
        if geneid == geneid_prev:
            if length >= count:
                count, protein = length, line[8][1:-1]
        else:
            for k, v in faa_dict.items():
                if protein in k:
                    result_sequence = k + v
                    outfile.write(result_sequence)
            count, protein = 0, line[8][1:-1]
        geneid_prev = geneid


if __name__ == '__main__':
    parser = ArgumentParser(description='FASTA filtering by ids from tabular file')
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str, help='FASTA file')
    group_required.add_argument('-t', '--tabular-file', type=str, help='Tabular (CSV) file from NCBI')
    group_required.add_argument('-o', '--output', type=str, help='outfile name')
    args = parser.parse_args()
    main()
