#!/usr/bin/env python3
__author__ = 'tomarovsky'

import csv
import argparse

input_file = 'all_mitos_annotations.1_based.seq'
output_file = 'all_mitos_annotations.1_based._seq'

def main():
    with open(args.input, 'r', newline='', encoding='utf-8') as infile, open(args.output, 'w', newline='', encoding='utf-8') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')

        for row in reader:
            if row and row[0].isdigit() and row[1].isdigit():
                # Increase the first and second columns by 1
                row[0] = str(int(row[0]) + 1)
                row[1] = str(int(row[1]) + 1)
            writer.writerow(row)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="script to convert mitos five-column (.seq) annotation coordinates from 0-based to 1-based")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str, default="all_mitos_annotations.seq", help="input .seq file")
    group_required.add_argument('-o', '--output', type=str, default="all_mitos_annotations.1_based.seq", help="output filename")
    args = parser.parse_args()
    main()



