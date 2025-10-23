#!/usr/bin/env python3
__author__ = 'tomarovsky'

import argparse
import gzip
import statistics
import sys
from collections import defaultdict


def open_file(filename):
    """Открывает обычный или .gz файл для чтения в текстовом режиме."""
    try:
        if filename.endswith(".gz"):
            return gzip.open(filename, 'rt', encoding='utf-8')
        else:
            return open(filename, 'r', encoding='utf-8')
    except FileNotFoundError:
        sys.stderr.write(f"Ошибка: Файл не найден {filename}\n")
        sys.exit(1)

def run_with_bed(vcf_file, bed_file, threshold):
    """Логика для режима VCF + BED."""

    valid_windows = defaultdict(list)
    kb_divisor = None

    # 1. Читаем BED и фильтруем окна
    with open_file(bed_file) as f_bed:
        for line in f_bed:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.split()
            if len(parts) < 4: continue # Простая проверка

            # Если 4-й столбец <= порога, окно подходит
            if float(parts[3]) <= threshold:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                # Сохраняем [start, end, snp_count]
                valid_windows[chrom].append([start, end, 0])

                # Определяем делитель по первому валидному окну
                if kb_divisor is None:
                    kb_divisor = (end - start) / 1000.0

    if kb_divisor is None or kb_divisor <= 0:
        sys.stderr.write("Ошибка: В BED не найдено валидных окон или размер окна 0.\n")
        sys.exit(1)

    # 2. Читаем VCF и считаем SNP
    total_snp_count = 0
    with open_file(vcf_file) as f_vcf:
        for line in f_vcf:
            if line.startswith("#"):
                continue

            total_snp_count += 1
            parts = line.split()
            chrom = parts[0]
            pos_0based = int(parts[1]) - 1 # VCF (1-based) -> BED (0-based)

            # Проверяем вхождение SNP в валидные окна для этой хромосомы
            if chrom in valid_windows:
                # window_data это [start, end, count]
                for window_data in valid_windows[chrom]:
                    if pos_0based >= window_data[0] and pos_0based < window_data[1]:
                        window_data[2] += 1 # Увеличиваем счетчик

    # 3. Собираем плотности
    snp_densities = []
    for chrom_list in valid_windows.values():
        for window_data in chrom_list:
            # window_data[2] - это snp_count
            snp_densities.append(window_data[2] / kb_divisor)

    return total_snp_count, snp_densities

def run_vcf_only(vcf_file):
    """Логика для режима "только VCF" (скользящие окна)."""

    window_size = 1000000
    step = 100000
    kb_divisor = window_size / 1000.0

    # { ('chr1', 0): 10, ('chr1', 100000): 12, ... }
    # Ключ = (хромосома, старт_окна)
    window_counts = defaultdict(int)
    total_snp_count = 0

    # 1. Читаем VCF и "раскидываем" SNP по окнам
    with open_file(vcf_file) as f_vcf:
        for line in f_vcf:
            if line.startswith("#"):
                continue

            total_snp_count += 1
            parts = line.split()
            chrom = parts[0]
            pos_0based = int(parts[1]) - 1

            # Определяем, в какие окна попадает SNP
            # Окно [s, s + 1M) содержит pos, если s <= pos < s + 1M
            # Это эквивалентно: pos - 1M < s <= pos

            # Самый правый (максимальный) старт окна, в который еще попадает SNP
            s_max = (pos_0based // step) * step

            # Самая левая граница для старта
            s_min_limit = pos_0based - window_size

            current_s = s_max
            while current_s > s_min_limit:
                key = (chrom, current_s)
                window_counts[key] += 1
                current_s -= step

    # 2. Собираем плотности
    snp_densities = [count / kb_divisor for count in window_counts.values()]

    return total_snp_count, snp_densities

def main():
    parser = argparse.ArgumentParser(description="Анализ плотности SNP.")
    parser.add_argument("-i", "--vcf", required=True, help="Путь к VCF файлу (может быть .gz).")
    parser.add_argument("-b", "--bed", required=False, default=None, help="(Опционально) Путь к BED файлу с окнами.")
    parser.add_argument("-f", "--filter_threshold", type=float, default=75000.0, help="Порог для 4-го столбца BED. (По умолчанию: 75000)")
    args = parser.parse_args()

    total_snps = 0
    densities = []

    try:
        if args.bed:
            # РЕЖИМ 1: VCF + BED
            total_snps, densities = run_with_bed(args.vcf, args.bed, args.filter_threshold)
        else:
            # РЕЖИМ 2: Только VCF
            total_snps, densities = run_vcf_only(args.vcf)

        # 4. Расчет финальной статистики и вывод
        if not densities:
            # Если VCF был пуст или BED не дал валидных окон
            mean_density = 0.0
            median_density = 0.0
        else:
            mean_density = statistics.mean(densities)
            median_density = statistics.median(densities)

        # Вывод в одну строку, разделенную табуляцией
        print(f"{total_snps}\t{mean_density:.4f}\t{median_density:.4f}")

    except Exception as e:
        sys.stderr.write(f"Критическая ошибка: {e}\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
