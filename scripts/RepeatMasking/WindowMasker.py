#!/usr/bin/env python3
__author__ = 'tomarovsky'
from Biocrutch.Routines.routine_functions import metaopen, metaoutput
import argparse
from sys import stdin


def main():
    if args.input == stdin:
        outprefix = 'windowmasker.asnb'
    else:
        if args.output:
            outprefix = args.output
        else:
            outprefix = ('.').join(args.input.split('.')[:-1])
    outprefix = metaoutput(outprefix, '.gff.gz')
    outfile = metaopen(outprefix, 'wt')
    print('INTERVAL to GFF converter...')
    with metaopen(args.input, 'r', buffering=args.buffering) as data:
        count = 0
        for l in data:
            if l.startswith('>'):
                line = l.strip().split('|')
                seq_name = line[1]
                # info = "_".join(line[2].split())
            else:
                line = l.strip().split(' ')
                start, stop = str(int(line[0]) + 1), line[2]
                gff_line = [
                    seq_name,
                    'WindowMasker',
                    'repeat',
                    start,
                    stop,
                    '.',
                    '.',
                    '.',
                    'ID=' + str(count) + '\n',
                ]  # INFO is just information from a FASTA file "+ ';INFO=' + info + '\n'"
                count += 1
                outfile.write('\t'.join(gff_line))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='windowmasker.asnb to GFF converter')
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str, default=stdin, help='input file or stdin')
    group_additional = parser.add_argument_group('Additional options')
    group_additional.add_argument('-b', '--buffering', type=int, default=None, help='Text buffering. Default = None')
    group_additional.add_argument('-o', '--output', type=str, help='output file prefix')
    args = parser.parse_args()
    main()
