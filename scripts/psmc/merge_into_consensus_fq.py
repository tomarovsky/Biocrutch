#!/usr/bin/env python3
__author__ = 'tomarovsky'

import sys
import textwrap

def fasta_parser(filename):
    header = None
    seq_parts = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith('>'):
                if header:
                    yield header, "".join(seq_parts)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
        if header:
            yield header, "".join(seq_parts)

def merge_to_fastq(fasta_file, qual_file):
    seq_gen = fasta_parser(fasta_file)
    qual_gen = fasta_parser(qual_file)

    wrapper = textwrap.TextWrapper(width=60, break_long_words=True, replace_whitespace=False)

    for (seq_header, seq), (qual_header, qual) in zip(seq_gen, qual_gen):
        if seq_header != qual_header:
            print(f"Critical error: Mismatched headers: {seq_header} != {qual_header}", file=sys.stderr)
            sys.exit(1)

        # Consensus header
        print(f"@{seq_header}")

        # Consensus sequence
        print("\n".join(wrapper.wrap(seq)))

        print("+")

        # Consensus quality
        print("\n".join(wrapper.wrap(qual)))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: merge_fastq.py <fasta_in> <qual_in>", file=sys.stderr)
        sys.exit(1)
    merge_to_fastq(sys.argv[1], sys.argv[2])
