#!/usr/bin/env python3
__author__ = 'tomarovsky'

import sys
import gzip

def split_consensus_fastq(fastq_gz_path, fasta_path, qual_path):
    state = None # 'seq' or 'qual'
    current_header = None
    current_seq = []
    current_qual = []

    try:
        if fastq_gz_path.endswith('.gz'):
             f_in = gzip.open(fastq_gz_path, 'rt')
        else:
            f_in = open(fastq_gz_path, 'r')

        with f_in, \
             gzip.open(fasta_path, 'wt') as f_out, \
             gzip.open(qual_path, 'wt') as q_out:

            def write_record():
                if current_header:
                    clean_header = current_header.split()[0]
                    f_out.write(f">{clean_header}\n")
                    f_out.write("".join(current_seq) + "\n")
                    q_out.write(f">{clean_header}\n")
                    q_out.write("".join(current_qual) + "\n")

            for line in f_in:
                line = line.strip()
                if not line:
                    continue

                if line.startswith('@'):
                    write_record()
                    current_header = line[1:]
                    current_seq = []
                    current_qual = []
                    state = 'seq'

                elif line == '+' and state == 'seq':
                    state = 'qual'

                elif state == 'seq':
                    current_seq.append(line)

                elif state == 'qual':
                    current_qual.append(line)

            write_record()

    except Exception as e:
        print(f"Error processing FASTQ: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    if len(sys.argv) != 4:
         print("Usage: split_consensus_fq.py <in_fastq> <out_fasta.gz> <out_qual.gz>", file=sys.stderr)
         sys.exit(1)
    split_consensus_fastq(sys.argv[1], sys.argv[2], sys.argv[3])
