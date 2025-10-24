#!/usr/bin/env python3
__author__ = 'tomarovsky'

import argparse
import gzip
import sys

import numpy as np
import pandas as pd

np_mean = np.mean
np_median = np.median


def open_file(filename):
    if filename.endswith(".gz"):
        return gzip.open(filename, 'rt', encoding='utf-8')
    else:
        return open(filename, 'r', encoding='utf-8')


def parse_scaffold_lengths_tsv(tsv_file):
    df = pd.read_csv(tsv_file, sep='\t', header=None,
                     names=['CHROM', 'LENGTH'], dtype={'CHROM': str, 'LENGTH': np.int64})
    return df.set_index('CHROM')['LENGTH'].to_dict()


def read_vcf_data(vcf_file, exclude_set):
    header_lines = 0
    with open_file(vcf_file) as f:
        for line in f:
            if line.startswith("##"):
                header_lines += 1
            elif line.startswith("#CHROM"):
                header_lines += 1
                break
            else:
                break

    df_vcf = pd.read_csv(
        vcf_file,
        sep='\t',
        skiprows=header_lines,
        header=None,
        usecols=[0, 1],
        names=['CHROM', 'POS'],
        dtype={'CHROM': str, 'POS': np.int32}
    )

    # to 0-based
    df_vcf['POS_0BASED'] = df_vcf['POS'] - 1

    # filter by exclude_set
    if exclude_set:
        df_vcf = df_vcf[~df_vcf['CHROM'].isin(exclude_set)]

    return df_vcf


def parse_bed_filter(bed_file, threshold):
    df_bed = pd.read_csv(
        bed_file,
        sep='\t',
        comment='#',
        header=None,
        usecols=[0, 1, 3],
        names=['CHROM', 'START', 'VALUE'],
        dtype={'CHROM': str, 'START': np.int32, 'VALUE': np.float64}
    )

    # filter by threshold
    df_filtered = df_bed[df_bed['VALUE'] <= threshold]

    # collect unique valid pairs (CHROM, START)
    valid_window_starts = set(zip(df_filtered['CHROM'], df_filtered['START']))

    return valid_window_starts


def run_window_analysis(df_vcf, scaffold_lengths, exclude_set, window_size, step, bed_filter_set=None):
    # Window size in kilobases (for density calculation)
    kb_divisor = window_size / 1000

    # 1. Generate all valid window starts
    window_data = []

    for chrom, length in scaffold_lengths.items():
        # Skip excluded scaffolds and those shorter than window size
        if chrom in exclude_set or length < window_size:
            continue

        # Generate all potential starts for sliding windows
        starts = np.arange(0, length - window_size + 1, step)

        for start in starts:
            key = (chrom, start)

            # Apply BED filter if provided
            if bed_filter_set is None or key in bed_filter_set:
                window_data.append((chrom, start, start + window_size))

    if not window_data:
        return []

    # Convert to DataFrame windows
    df_windows = pd.DataFrame(window_data, columns=['CHROM', 'START', 'END'])

    # Initialize counters
    # Important: use df_windows['CHROM'].astype('category') for faster grouping
    df_windows['CHROM'] = df_windows['CHROM'].astype('category')
    window_counts = pd.Series(0, index=df_windows.index, dtype=np.int32)
    grouped_vcf = df_vcf.groupby('CHROM')

    # 2. Count SNPs in windows
    for chrom in df_windows['CHROM'].unique():
        try:
            snp_subset = grouped_vcf.get_group(chrom)
        except KeyError:
            continue

        # Select windows for current chromosome
        window_subset = df_windows[df_windows['CHROM'] == chrom]

        # Convert SNP positions to NumPy array for efficient comparison
        snp_positions = snp_subset['POS_0BASED'].values

        # Create 2D arrays for vectorized comparison
        # Create arrays of window starts and ends for current chromosome
        starts = window_subset['START'].values[:, None]  # Column vector
        ends = window_subset['END'].values[:, None]      # Column vector

        # Count how many SNPs fall into each window.
        # (N_windows, N_snps) array of boolean values.
        # Boolean mask: (start <= POS) & (POS < end)
        mask = (snp_positions >= starts) & (snp_positions < ends)

        # Sum along rows (window), to get the number of SNPs in each window
        counts = mask.sum(axis=1)

        # Update Series counters.
        window_counts.loc[window_subset.index] = counts

    # 3. Calculate density (SNPs/kb)
    snp_densities = (window_counts / kb_divisor).tolist()

    return snp_densities


def main():
    parser = argparse.ArgumentParser(description="SNP density analysis using Pandas for performance.")
    parser.add_argument("-i", "--vcf", required=True, help="Path to the VCF file.")
    parser.add_argument("-l", "--scaffold_lengths", required=True,
                        help="Path to the TSV file with scaffold lengths (scaffold_name\\tlength).")
    parser.add_argument("-b", "--bed", required=False, default=None,
                        help="(Optional) Path to a BED file. If specified, it is used to filter sliding windows based on the 4th column.")
    parser.add_argument("-t", "--threshold", type=float, default=75000,
                        help="Threshold for the 4th BED column (only used when --bed is specified).")
    parser.add_argument("-w", "--window_size", type=int, default=1000000,
                        help="Sliding window size in bp. (Default: 1000000)")
    parser.add_argument("-s", "--step_size", type=int, default=100000,
                        help="Sliding window step size in bp. (Default: 100000)")
    parser.add_argument("-e", "--exclude", required=False, nargs='+', help="List of scaffolds/chromosomes to exclude.")

    args = parser.parse_args()
    exclude_set = set(args.exclude) if args.exclude else set()

    try:
        scaffold_lengths = parse_scaffold_lengths_tsv(args.scaffold_lengths)
        df_vcf = read_vcf_data(args.vcf, exclude_set)
        total_snps_final = len(df_vcf)

        bed_filter_set = None
        if args.bed:
            bed_filter_set = parse_bed_filter(args.bed, args.threshold)

            if not bed_filter_set:
                sys.stderr.write("Warning: BED file specified, but no windows passed the filter. Densities will be 0.0.\n")

        densities = run_window_analysis(
            df_vcf,
            scaffold_lengths,
            exclude_set,
            args.window_size,
            args.step_size,
            bed_filter_set
        )

        if densities:
            mean_density = np_mean(densities)
            median_density = np_median(densities)

        # Total_SNPs Mean_Density Median_Density
        print(f"{total_snps_final}\t{mean_density:.4f}\t{median_density:.4f}")

    except Exception as e:
        sys.stderr.write(f"Critical error: {e}\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
