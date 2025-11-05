#!/usr/bin/env python3
import argparse
import gzip
import pandas as pd
import numpy as np

def load_bed(bed_file):
    """Load BED file into a pandas DataFrame and group by scaffold."""
    bed = pd.read_csv(bed_file, sep='\t', header=None, usecols=[0, 1, 2],
                     names=['scaffold', 'start', 'end'])
    return {scaffold: list(zip(group['start'], group['end']))
            for scaffold, group in bed.groupby('scaffold')}

def parse_consensus_fastq(handle):
    """
    Generator: yields (header, seq_lines, qual_lines)
    Assumes that for each scaffold:
        header line -> sequence lines -> '+' line -> quality lines
        number of sequence lines == number of quality lines
    """
    while True:
        header_line = handle.readline()
        if not header_line:
            break
        header = header_line.rstrip()[1:]  # remove '@'

        # read sequence lines until '+'
        seq_lines = []
        while True:
            line = handle.readline()
            if not line:
                raise ValueError(f"Unexpected end of file while reading sequence for {header}")
            line = line.rstrip()
            if line.startswith('+'):
                break
            seq_lines.append(line)

        # read same number of quality lines as sequence lines
        qual_lines = []
        for _ in seq_lines:
            qline = handle.readline()
            if not qline:
                raise ValueError(f"Unexpected end of file while reading quality for {header}")
            qual_lines.append(qline.rstrip())

        yield header, ''.join(seq_lines), ''.join(qual_lines)

def mask_with_numpy(seq, qual, regions):
    """Fast masking of consensus FASTQ using NumPy"""
    seq_arr = np.frombuffer(seq.encode('ascii'), dtype='S1').copy()
    qual_arr = np.frombuffer(qual.encode('ascii'), dtype='S1').copy()

    for start, end in regions:
        seq_arr[start:end] = b'n'
        qual_arr[start:end] = b'!'

    return seq_arr.tobytes().decode('ascii'), qual_arr.tobytes().decode('ascii')

def main():
    parser = argparse.ArgumentParser(description="Masking of consensus FASTQ")
    parser.add_argument("-i", "--fastq", required=True, help="Input consensus FQ.gz")
    parser.add_argument("-m", "--bed", required=True, help="BED file with regions to mask")
    parser.add_argument("-o", "--output", required=True, help="Output masked consensus FQ.gz")
    args = parser.parse_args()

    bed_dict = load_bed(args.bed)

    open_func = gzip.open if args.fastq.endswith(".gz") else open
    out_func = gzip.open if args.output.endswith(".gz") else open

    with open_func(args.fastq, "rt") as fin, out_func(args.output, "wt") as fout:
        for header, seq, qual in parse_consensus_fastq(fin):
            if header in bed_dict:
                seq, qual = mask_with_numpy(seq, qual, bed_dict[header])
            fout.write(f"@{header}\n")
            for i in range(0, len(seq), 60):
                fout.write(seq[i:i+60] + '\n')
            fout.write("+\n")
            for i in range(0, len(qual), 60):
                fout.write(qual[i:i+60] + '\n')

if __name__ == "__main__":
    main()
