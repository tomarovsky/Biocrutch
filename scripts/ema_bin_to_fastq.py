from Biocrutch.Routines.routine_functions import metaopen
import argparse
from sys import stdin


def main():
    if args.input == '-':
        with metaopen(args.prefix + "_barcodes.txt.gz", 'wt') as barcode_fd, metaopen(args.prefix + "_ema-bin-all_1.fastq.gz", 'wt') as forward_fd, metaopen(args.prefix + "_ema-bin-all_2.fastq.gz", 'wt') as reverse_fd:
            for line in stdin:
                line = line.split()
                barcode_fd.write(line[0] + '\n')
                forward_fd.write(line[1] + '\n' + line[2] + '\n' + '+\n' + line[3] + '\n')
                reverse_fd.write(line[1] + '\n' + line[4] + '\n' + '+\n' + line[5] + '\n')
    else:
        with metaopen(args.input, "r", args.buffering) as fd, metaopen(args.prefix + "_barcodes.txt.gz", 'wt') as barcode_fd, metaopen(args.prefix + "_ema-bin-all_1.fastq.gz", 'wt') as forward_fd, metaopen(args.prefix + "_ema-bin-all_2.fastq.gz", 'wt') as reverse_fd:
            for line in fd:
                line = line.split()
                barcode_fd.write(line[0] + '\n')
                forward_fd.write(line[1] + '\n' + line[2] + '\n' + '+\n' + line[3] + '\n')
                reverse_fd.write(line[1] + '\n' + line[4] + '\n' + '+\n' + line[5] + '\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ema-bin-xxx to forward_1.fastq, reverse_2.fastq and barcodes.txt")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str,
                             help="input file ema-bin-xxx (use '-' for STDIN)")
    group_required.add_argument('-p', '--prefix', type=str,
                             help="output files prefix")
    group_additional = parser.add_argument_group('Additional options')
    group_additional.add_argument('-b', '--buffering',
                             metavar='INT', type=int, 
                             default=None, help="Text buffering. Default = None")
    args = parser.parse_args()
    main()


