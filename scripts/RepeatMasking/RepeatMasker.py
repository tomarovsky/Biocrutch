#!/usr/bin/env python3
__author__ = 'tomarovsky'
from Biocrutch.Routines.routine_functions import metaopen
import argparse
from sys import stdin

def main():
    with metaopen(args.input, "r", buffering=args.buffering) as data:
        for line in data:
            print(line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="repeatmasker.out to GFF converter")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str, default=stdin,
                                help="input file or stdin")
    group_additional = parser.add_argument_group('Additional options')
    group_additional.add_argument('-b', '--buffering', type=int, 
                                  default=None, help="Text buffering. Default = None")
    args = parser.parse_args()
    main()