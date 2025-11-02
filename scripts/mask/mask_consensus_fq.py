#!/usr/bin/env python3
import sys
import gzip

def read_bed(bed_file):
    """Read a BED file and return a dictionary of scaffold -> list of (start, end) regions."""
    mask_regions = {}
    with open(bed_file) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, *rest = line.strip().split()
            start, end = int(start), int(end)
            mask_regions.setdefault(chrom, []).append((start, end))
    return mask_regions


def mask_sequence(seq, qual, regions):
    """Mask a sequence and quality scores in specified regions."""
    if not regions:
        return seq, qual
    seq = list(seq)
    qual = list(qual)
    for start, end in regions:
        # BED coordinates are 0-based, end is not inclusive
        for i in range(start, min(end, len(seq))):
            seq[i] = 'n'
            qual[i] = '!'
    return ''.join(seq), ''.join(qual)


def read_consensus_fastq(fastq_path):
    """Generator for reading consensus FASTQ (one scaffold at a time)."""
    open_func = gzip.open if fastq_path.endswith(".gz") else open
    with open_func(fastq_path, "rt") as f:
        while True:
            header = f.readline()
            if not header:
                break
            if not header.startswith('@'):
                raise ValueError(f"Unexpected FASTQ header line: {header.strip()}")
            name = header[1:].strip()

            # читаем последовательность до '+'
            seq_lines = []
            for line in f:
                if line.startswith('+'):
                    break
                seq_lines.append(line.strip())
            seq = ''.join(seq_lines)

            # read quality
            qual_lines = []
            while len(''.join(qual_lines)) < len(seq):
                line = f.readline()
                if not line:
                    break
                qual_lines.append(line.strip())
            qual = ''.join(qual_lines)

            yield name, seq, qual


def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} CONSENSUS_FASTQ(.gz) MASK_BED OUTPUT_FASTQ(.gz)", file=sys.stderr)
        sys.exit(1)

    fastq_path, bed_path, out_path = sys.argv[1:]

    mask_regions = read_bed(bed_path)

    open_func_out = gzip.open if out_path.endswith(".gz") else open
    with open_func_out(out_path, "wt") as out_f:
        for name, seq, qual in read_consensus_fastq(fastq_path):
            regions = mask_regions.get(name, [])
            masked_seq, masked_qual = mask_sequence(seq, qual, regions)

            # write masked sequence and quality
            out_f.write(f"@{name}\n")
            for i in range(0, len(masked_seq), 60):
                out_f.write(masked_seq[i:i+60] + "\n")
            out_f.write("+\n")
            for i in range(0, len(masked_qual), 60):
                out_f.write(masked_qual[i:i+60] + "\n")


if __name__ == "__main__":
    main()
