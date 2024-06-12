#!/usr/bin/env python3
__author__ = 'tomarovsky'

import argparse
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq

def parse_exclude_list(exclude_file):
    """
    Parses the file containing the lists of scaffolds to exclude, coordinates for trimming,
    and duplicated scaffolds.
    :param exclude_file: Path to the text file with sections "Exclude:", "Trim:", and "Duplicated:".
    :return: A set of scaffolds to exclude, a dictionary with coordinates for trimming, and a set of duplicated scaffolds.
    """
    exclude_scaffolds = set()
    trim_info = {}
    duplicated_scaffolds = set()

    if exclude_file:
        with open(exclude_file, "r") as file:
            exclude_section = False
            trim_section = False
            duplicated_section = False
            for line in file:
                if line.startswith("Exclude:"):
                    exclude_section = True
                    trim_section = False
                    duplicated_section = False
                    continue
                if line.startswith("Trim:"):
                    trim_section = True
                    exclude_section = False
                    duplicated_section = False
                    continue
                if line.startswith("Duplicated:"):
                    duplicated_section = True
                    exclude_section = False
                    trim_section = False
                    continue
                if exclude_section:
                    if not line.strip() or line.startswith("Sequence"):
                        continue
                    scaffold_name = line.split()[0]
                    exclude_scaffolds.add(scaffold_name)
                if trim_section:
                    if not line.strip() or line.startswith("Sequence"):
                        continue
                    parts = line.split()
                    # if len(parts) < 4:
                    #     continue
                    scaffold_name = parts[0]
                    coordinates = parts[2]
                    start, stop = map(int, coordinates.split(".."))
                    if scaffold_name not in trim_info:
                        trim_info[scaffold_name] = []
                    trim_info[scaffold_name].append((start, stop))
                if duplicated_section:
                    if not line.strip() or line.startswith("#") or line.startswith("Sequence"):
                        continue
                    scaffolds = line.split()[:-2]
                    for scaffold in scaffolds[1:]:
                        if scaffold.startswith("RC("):
                            scaffold = scaffold[3:-1]
                        duplicated_scaffolds.add(scaffold)
    return exclude_scaffolds, trim_info, duplicated_scaffolds

def trim_sequences(sequence, trim_coordinates):
    """
    Replaces sequence regions with "N" based on the given coordinates.
    :param sequence: SeqRecord object to process.
    :param trim_coordinates: List of tuples with coordinates (start, stop).
    :return: Processed SeqRecord object.
    """
    mutable_seq = MutableSeq(str(sequence.seq))
    for start, stop in trim_coordinates:
        mutable_seq[start-1:stop] = 'N' * (stop - start + 1)
    sequence.seq = Seq(str(mutable_seq))
    return sequence

def filter_and_trim_contigs(input_fasta, output_fasta, min_length, exclude_scaffolds, trim_info, duplicated_scaffolds):
    """
    Filters out contigs that are <= min_length, listed in the exclusion set,
    replaces regions with "N" based on the coordinates, and removes duplicated scaffolds.
    :param input_fasta: Path to the input FASTA file.
    :param output_fasta: Path to the output FASTA file.
    :param min_length: Minimum contig length to keep.
    :param exclude_scaffolds: Set of scaffolds to exclude.
    :param trim_info: Dictionary with coordinates for trimming.
    :param duplicated_scaffolds: Set of duplicated scaffolds to exclude.
    """
    with open(input_fasta, "r") as input_handle, open(output_fasta, "w") as output_handle:
        sequences = SeqIO.parse(input_handle, "fasta")
        filtered_sequences = (seq for seq in sequences if len(seq) > min_length and seq.id not in exclude_scaffolds and seq.id not in duplicated_scaffolds)
        for seq in filtered_sequences:
            if seq.id in trim_info:
                seq = trim_sequences(seq, trim_info[seq.id])
            SeqIO.write(seq, output_handle, "fasta")

def main():
    parser = argparse.ArgumentParser(description="NCBI contamination filtration.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file")
    parser.add_argument("-m", "--min_length", type=int, default=200, help="Minimal contig length")
    parser.add_argument("-c", "--contamination", help="NCBI contamination file", default=None)

    args = parser.parse_args()

    exclude_scaffolds, trim_info, duplicated_scaffolds = parse_exclude_list(args.contamination)
    filter_and_trim_contigs(args.input, args.output, args.min_length, exclude_scaffolds, trim_info, duplicated_scaffolds)

if __name__ == "__main__":
    main()

