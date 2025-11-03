#!/usr/bin/env python3
import argparse
import gzip
import pandas as pd
import numpy as np

def wrap_string(s, width=60):
    """Wraps a string into lines of specified width."""
    return '\n'.join(s[i:i+width] for i in range(0, len(s), width))

def load_bed(bed_file):
    """Load BED file into a pandas DataFrame and group by scaffold."""
    bed = pd.read_csv(bed_file, sep='\t', header=None, usecols=[0, 1, 2],
                     names=['scaffold', 'start', 'end'])
    return {scaffold: list(zip(group['start'], group['end']))
            for scaffold, group in bed.groupby('scaffold')}

def parse_consensus_fastq(handle):
    """Generator: yields (header, seq, qual)."""
    header, seq_lines, qual_lines = None, [], []
    reading_seq, reading_qual = False, False

    for line in handle:
        line = line.rstrip('\n')
        if line.startswith('@'):
            if header:
                yield header, ''.join(seq_lines), ''.join(qual_lines)
            header = line[1:].strip()
            seq_lines, qual_lines = [], []
            reading_seq, reading_qual = True, False
        elif line.startswith('+'):
            reading_seq, reading_qual = False, True
        else:
            if reading_seq:
                seq_lines.append(line)
            elif reading_qual:
                qual_lines.append(line)

    if header:
        yield header, ''.join(seq_lines), ''.join(qual_lines)

def mask_with_numpy(seq, qual, regions):
    """Fast masking of consensus FASTQ using NumPy"""
    seq_ba = bytearray(seq.encode('ascii'))
    qual_ba = bytearray(qual.encode('ascii'))

    seq_arr = np.frombuffer(seq_ba, dtype='S1')
    qual_arr = np.frombuffer(qual_ba, dtype='S1')

    mask = np.zeros(len(seq_arr), dtype=bool)
    for start, end in regions:
        mask[start:end] = True

    seq_arr[mask] = b'n'
    qual_arr[mask] = b'!'

    return seq_ba.decode('ascii'), qual_ba.decode('ascii')


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
            fout.write(wrap_string(seq, 60) + '\n')
            fout.write("+\n")
            fout.write(wrap_string(qual, 60) + '\n')

if __name__ == "__main__":
    main()
