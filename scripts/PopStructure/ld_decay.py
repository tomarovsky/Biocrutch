#!/usr/bin/env python3

import gzip
import sys
import argparse
import numpy as np
from collections import defaultdict


def getOptionValue(option):
    if option in sys.argv:
        optionPos = sys.argv.index(option)
        optionValue = sys.argv[optionPos + 1]
        return optionValue
    else:
        print('\nWarning, option', option, 'not specified.\n', file=sys.stderr)


parser = argparse.ArgumentParser(description='Calculate LD decay from gzipped LD file.')
parser.add_argument('-i', '--input', type=str, help='Input gzipped file name', required=True)
parser.add_argument('-o', '--output', type=str, help='Output prefix', required=True)
parser.add_argument('-b', '--bins', type=int, default=100, help='Number of bins for averaging')

args = parser.parse_args()

fileName = args.input
prefix = args.output
num_bins = args.bins

f = gzip.open(fileName, 'rt')

# skip the header line
next(f)

chr_distance_ld = defaultdict(lambda: defaultdict(list))

for i, line in enumerate(f):
    line = line.strip().split()
    distance = abs(int(line[1]) - int(line[4]))
    chr = line[0]

    chr_distance_ld[chr][distance].append(float(line[6]))

    if (i + 1) % 100000000 == 0:
        print(f'Read {i + 1} entries.', file=sys.stderr)

print('Calculating average distances.', file=sys.stderr)

decay_output = open(f'{prefix}.ld_decay', 'w')
decay_output.write('chr\tdistance\tavg_R2\tstd\n')
decay_output_bins = open(f'{prefix}.ld_decay_bins', 'w')
decay_output_bins.write('chr\tdistance\tavg_R2\tstd\n')

for chr, distances in chr_distance_ld.items():
    distance_bin = defaultdict(list)
    for distance, lds in distances.items():
        mean = np.mean(lds)
        std = np.std(lds)
        decay_output.write(f'{chr}\t{distance}\t{mean:.6f}\t{std:.6f}\n')
        bin_size = int(round(distance / num_bins, 0)) * num_bins
        distance_bin[bin_size].extend(lds)
    for bin, lds in distance_bin.items():
        mean = np.mean(lds)
        std = np.std(lds)
        decay_output_bins.write(f'{chr}\t{bin}\t{mean:.6f}\t{std:.6f}\n')

decay_output_bins.close()
decay_output.close()

f.close()
