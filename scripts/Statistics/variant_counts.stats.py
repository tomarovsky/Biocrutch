#!/usr/bin/env python3
__author__ = 'tomarovsky'

import argparse
import gzip
import statistics
import sys
from collections import defaultdict


def open_file(filename):
    if filename.endswith(".gz"):
        return gzip.open(filename, 'rt', encoding='utf-8')
    else:
        return open(filename, 'r', encoding='utf-8')


def parse_scaffold_lengths_tsv(tsv_file):
    scaffold_lengths = {}
    with open_file(tsv_file) as f:
        for line in f:
            parts = line.split()
            chrom, length = parts[0], int(parts[1])
            scaffold_lengths[chrom] = length
    return scaffold_lengths


def run_with_bed(vcf_file, bed_file, threshold, exclude_set):
    """Логика для режима VCF + BED."""
    valid_windows = defaultdict(list)
    kb_divisor = None

    # 1. Читаем BED, фильтруем окна по --exclude и --filter_threshold
    with open_file(bed_file) as f_bed:
        for line in f_bed:
            parts = line.split()

            chrom = parts[0]
            if chrom in exclude_set or float(parts[3]) > threshold:
                continue

            start, end = int(parts[1]), int(parts[2])
            valid_windows[chrom].append([start, end, 0])

            if kb_divisor is None:
                kb_divisor = (end - start) / 1000.0

    # 2. Читаем VCF и считаем SNP
    total_snp_count = 0
    with open_file(vcf_file) as f_vcf:
        for line in f_vcf:
            if line.startswith("#"):
                continue

            parts = line.split()
            chrom = parts[0]

            if chrom in exclude_set:
                continue

            total_snp_count += 1
            pos_0based = int(parts[1]) - 1

            if chrom in valid_windows:
                for window_data in valid_windows[chrom]:
                    # Интервал [start, end)
                    if pos_0based >= window_data[0] and pos_0based < window_data[1]:
                        window_data[2] += 1

    # 3. Собираем плотности (SNP/kb)
    snp_densities = [window[2] / kb_divisor for chrom_list in valid_windows.values() for window in chrom_list]

    return total_snp_count, snp_densities


def run_vcf_only(vcf_file, scaffold_lengths, exclude_set):
    """Логика для режима VCF-only (скользящие окна 1Мб/100кб) с учетом длин скаффолдов."""
    window_size = 1000000
    step = 100000
    kb_divisor = window_size / 1000

    # 1. Генерируем список ВАЛИДНЫХ окон
    window_counts = {}
    valid_scaffolds_for_windowing = set()

    for chrom, length in scaffold_lengths.items():
        if chrom in exclude_set or length < window_size:
            continue

        valid_scaffolds_for_windowing.add(chrom)

        # Создаем окна, которые ПОЛНОСТЬЮ помещаются на скаффолде.
        current_start = 0
        while (current_start + window_size) <= length:
            window_counts[(chrom, current_start)] = 0
            current_start += step

    # 2. Читаем VCF, считаем Total_SNPs и распределяем по валидным окнам
    total_snp_count = 0
    with open_file(vcf_file) as f_vcf:
        for line in f_vcf:
            if line.startswith("#"):
                continue

            parts = line.split()
            chrom = parts[0]

            if chrom in exclude_set:
                continue

            total_snp_count += 1

            if chrom not in valid_scaffolds_for_windowing:
                continue

            pos_0based = int(parts[1]) - 1

            # Определяем, в какие перекрывающиеся окна попадает SNP
            s_max = (pos_0based // step) * step
            s_min_limit = pos_0based - window_size + 1

            current_s = s_max
            while current_s >= 0 and current_s >= s_min_limit:
                key = (chrom, current_s)
                # Увеличиваем счетчик, только если окно является "валидным" (полностью влезает в скаффолд)
                if key in window_counts:
                    window_counts[key] += 1
                current_s -= step

    # 3. Собираем плотности (SNP/kb)
    snp_densities = [count / kb_divisor for count in window_counts.values()]

    return total_snp_count, snp_densities


def main():
    parser = argparse.ArgumentParser(description="Анализ плотности SNP.")
    parser.add_argument("-i", "--vcf", required=True, help="path to VCF file")
    parser.add_argument("-l", "--scaffold_length_file", required=True, help=".len file")
    parser.add_argument("-b", "--bed", required=False, default=None, help="path to BED counts file")
    parser.add_argument("-t", "--filter_threshold", type=float, default=None, help="threshold for filtering")
    parser.add_argument("--exclude", required=False, nargs='+', help="list of scaffolds/chromosomes to exclude")

    args = parser.parse_args()
    exclude_set = set(args.exclude) if args.exclude else set()

    total_snps = 0
    densities = []

    try:
        # Считывание длин скаффолдов обязательно
        scaffold_lengths = parse_scaffold_lengths_tsv(args.scaffold_length_file)

        if args.bed:
            total_snps, densities = run_with_bed(
                args.vcf, args.bed, args.filter_threshold, exclude_set
            )
        else:
            total_snps, densities = run_vcf_only(args.vcf, scaffold_lengths, exclude_set)

        if densities:
            mean_density = statistics.mean(densities)
            median_density = statistics.median(densities)
        else:
            mean_density = 0.0
            median_density = 0.0

        # Формат вывода: Total_SNPs Mean_Density Median_Density
        print(f"{total_snps}\t{mean_density:.4f}\t{median_density:.4f}")

    except Exception as e:
        sys.stderr.write(f"Критическая ошибка: {e}\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
