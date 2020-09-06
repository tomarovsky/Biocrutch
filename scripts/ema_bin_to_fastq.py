from Biocrutch.Routines.routine_functions import metaopen
import argparse
from sys import stdin


def main():
    barcode_fd = metaopen(args.prefix + "_barcodes.txt.gz", 'wt')
    forward_fd = metaopen(args.prefix + "_ema-bin-all_1.fastq.gz", 'wt')
    reverse_fd = metaopen(args.prefix + "_ema-bin-all_2.fastq.gz", 'wt')
    if args.input == '-':
        data = stdin
    elif args.compressed: #for file.gz
        data = metaopen(args.input, "rt", args.buffering)
    else:
        data = metaopen(args.input, "r", args.buffering)
    for line in data:
        line = line.split()
        barcode_fd.write(line[0] + '\n')
        forward_fd.write(line[1] + '\n' + line[2] + '\n' + '+\n' + line[3] + '\n')
        reverse_fd.write(line[1] + '\n' + line[4] + '\n' + '+\n' + line[5] + '\n')
    barcode_fd.close()
    forward_fd.close()
    reverse_fd.close()
    if args.input != '-':
        data.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ema-bin-xxx to forward_1.fastq.gz, reverse_2.fastq.gz and barcodes.txt.gz")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str,
                             help="input file ema-bin-xxx (use '-' for STDIN)")
    group_required.add_argument('-p', '--prefix', type=str,
                             help="output files prefix")
    group_additional = parser.add_argument_group('Additional options')
    group_additional.add_argument('-c', '--compressed', action="store_true",
                             help="Flag. if your input is compressed (file.gz)")
    group_additional.add_argument('-b', '--buffering',
                             metavar='INT', type=int, 
                             default=None, help="Text buffering. Default = None")
    args = parser.parse_args()
    main()


