#!/usr/bin/env python3
__author__ = 'tomarovsky'

import sys
import textwrap
import gzip

def fasta_parser(filename):
    header = None
    seq_parts = []
    try:
        if filename.endswith('.gz'):
            f = gzip.open(filename, 'rt')
        else:
            f = open(filename, 'r')

        with f:
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
    except Exception as e:
        print(f"Parsing error in {filename}: {e}", file=sys.stderr)
        sys.exit(1)


def merge_to_fastq(fasta_file, qual_file, output_fastq_gz):
    seq_gen = fasta_parser(fasta_file)
    qual_gen = fasta_parser(qual_file)

    wrapper = textwrap.TextWrapper(width=60, break_long_words=True, replace_whitespace=False)

    try:
        with gzip.open(output_fastq_gz, 'wt') as f_out:
            for (seq_header, seq), (qual_header, qual) in zip(seq_gen, qual_gen):
                if seq_header != qual_header:
                    print(f"Critical error: Mismatched headers: {seq_header} != {qual_header}", file=sys.stderr)
                    sys.exit(1)

                f_out.write(f"@{seq_header}\n")
                f_out.write("\n".join(wrapper.wrap(seq)) + "\n")
                f_out.write("+\n")
                f_out.write("\n".join(wrapper.wrap(qual)) + "\n")

    except Exception as e:
        print(f"Writing error in {output_fastq_gz}: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: merge_into_consensus_fq.py <fasta_in.gz> <qual_in.gz> <output_fastq_gz>", file=sys.stderr)
        sys.exit(1)
    merge_to_fastq(sys.argv[1], sys.argv[2], sys.argv[3])
