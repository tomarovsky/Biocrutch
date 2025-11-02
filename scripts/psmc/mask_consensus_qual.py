#!/usr/bin/env python3
__author__ = 'tomarovsky'

import sys

def fasta_parser(filename):
    header = None
    seq_parts = []
    try:
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
    except Exception as e:
        print(f"An error occurred while parsing {filename}: {e}", file=sys.stderr)
        sys.exit(1)

def mask_quality(masked_fasta_file, original_qual_file, output_qual_file):
    seq_gen = fasta_parser(masked_fasta_file)
    qual_gen = fasta_parser(original_qual_file)

    with open(output_qual_file, 'w') as q_out:
        for (seq_header, masked_seq), (qual_header, original_qual) in zip(seq_gen, qual_gen):
            if seq_header != qual_header:
                print(f"Critical error: Mismatched headers in masked FASTA and original QUAL: {seq_header} != {qual_header}", file=sys.stderr)
                sys.exit(1)

            if len(masked_seq) != len(original_qual):
                print(f"Critical error: Sequence/Quality length mismatch for {seq_header}. Seq length: {len(masked_seq)}, Qual length: {len(original_qual)}", file=sys.stderr)
                sys.exit(1)

            masked_qual_list = list(original_qual)

            for i in range(len(masked_seq)):
                if masked_seq[i].lower() == 'n':
                    masked_qual_list[i] = '!'

            masked_qual = "".join(masked_qual_list)

            # Write the masked quality record in FASTA format
            q_out.write(f">{qual_header}\n")
            q_out.write(masked_qual + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: mask_qual.py <masked_fasta> <original_qual> <output_masked_qual>", file=sys.stderr)
        sys.exit(1)
    mask_quality(sys.argv[1], sys.argv[2], sys.argv[3])
