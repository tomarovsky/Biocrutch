from Biocrutch.Routines.routine_functions import metaopen
import argparse
import os


def main():
    for file in args.input:
        fh = metaopen(file, "r", args.buffering)
        for line in fh:
            line = line.split()
            with open("barcodes.txt", 'a') as f:
                f.write(line[0] + '\n')
            with open(os.path.splitext(file)[0] + "_1.fastq", 'a') as f:
                f.write(line[1] + '\n' + line[2] + '\n' + '+\n' + line[3] + '\n')
            with open(os.path.splitext(file)[0] + "_2.fastq", 'a') as f:
                f.write(line[1] + '\n' + line[4] + '\n' + '+\n' + line[5] + '\n')
        fh.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ema-bin-xxx to forward_1.fastq, reverse_2.fastq and barcodes.txt")
    group_required = parser.add_argument_group('Options')
    group_required.add_argument('-i', '--input', type=str,
                            nargs='+', help="input file ema-bin-xxx")
    group_required.add_argument('-b', '--buffering',
                             metavar='INT', type=int, 
                             default=None, help="Text buffering. Default = None")
    args = parser.parse_args()
    main()


