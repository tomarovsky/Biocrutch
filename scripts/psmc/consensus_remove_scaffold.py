#!/usr/bin/env python3
import gzip
import argparse

def remove_fq_block_by_id(input_file, output_file, remove_id):
    open_func = gzip.open if input_file.endswith(".gz") else open
    with open_func(input_file, 'rt') as f, open_func(output_file, 'wt') as out:
        while True:
            header = f.readline()
            if not header:
                break
            header = header.rstrip()
            # read sequence
            seq_lines = []
            while True:
                line = f.readline().rstrip()
                if line == '+':
                    plus_line = line
                    break
                seq_lines.append(line)
            # read quality — same number of lines as seq_lines
            qual_lines = [f.readline().rstrip() for _ in range(len(seq_lines))]

            if header[1:].rstrip() == remove_id:
                continue  # skip block
            # write block
            out.write(header + "\n")
            out.write("\n".join(seq_lines) + "\n")
            out.write(plus_line + "\n")
            out.write("\n".join(qual_lines) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove scaffold from consensus FQ by ID")
    parser.add_argument("-i", "--input", required=True, help="Input FQ(.gz) file")
    parser.add_argument("-o", "--output", required=True, help="Output FQ(.gz) file")
    parser.add_argument("-s", "--remove_id", required=True, help="ID to remove (without @)")
    args = parser.parse_args()

    remove_fq_block_by_id(args.input, args.output, args.remove_id)
