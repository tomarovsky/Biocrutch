#!/usr/bin/env python3
__author__ = 'tomarovsky'

import argparse
import gzip
import statistics
import sys

import numpy as np
import pandas as pd


def open_file(filename):
    """Возвращает файловый дескриптор или gzip.open."""
    if filename.endswith(".gz"):
        return gzip.open(filename, 'rt', encoding='utf-8')
    else:
        return open(filename, 'r', encoding='utf-8')


def parse_scaffold_lengths_tsv(tsv_file):
    """Парсит TSV файл с длинами скаффолдов с помощью Pandas."""
    try:
        # Читаем TSV
        df = pd.read_csv(tsv_file, sep='\s+', header=None, names=['CHROM', 'LENGTH'], dtype={'CHROM': str, 'LENGTH': np.int64})
        if df.empty:
            raise ValueError(f"Файл {tsv_file} не содержит валидных данных.")
        return df.set_index('CHROM')['LENGTH'].to_dict()
    except Exception as e:
        sys.stderr.write(f"Ошибка при чтении TSV с длинами: {e}\n")
        sys.exit(1)


def read_vcf_data(vcf_file, exclude_set):
    """
    Читает VCF файл (только данные, без заголовка) в Pandas DataFrame.
    Сразу фильтрует по exclude_set.
    """
    header_lines = 0
    # Находим конец заголовка VCF
    with open_file(vcf_file) as f:
        for line in f:
            if line.startswith("#"):
                header_lines += 1
            else:
                break

    try:
        df_vcf = pd.read_csv(
            vcf_file,
            sep='\t',
            skiprows=header_lines - 1,
            header=0,
            usecols=[0, 1],
            names=['CHROM', 'POS'],
            dtype={'CHROM': str, 'POS': np.int32}
        )
        # 0-based позиция
        df_vcf['POS_0BASED'] = df_vcf['POS'] - 1

        # Фильтрация по exclude_set
        if exclude_set:
             df_vcf = df_vcf[~df_vcf['CHROM'].isin(exclude_set)]

        return df_vcf
    except Exception as e:
        sys.stderr.write(f"Ошибка при чтении VCF данных: {e}\n")
        sys.exit(1)


def parse_bed_filter(bed_file, threshold):
    """Читает BED и возвращает множество валидных (CHROM, START) окон, прошедших фильтр по 4-му столбцу."""
    valid_window_starts = set()
    try:
        df_bed = pd.read_csv(
            bed_file,
            sep='\t',
            comment='#',
            header=None,
            usecols=[0, 1, 3],
            names=['CHROM', 'START', 'VALUE'],
            dtype={'CHROM': str, 'START': np.int32, 'VALUE': np.float64}
        )

        # Фильтруем по порогу VALUE <= threshold
        df_filtered = df_bed[df_bed['VALUE'] <= threshold]

        # Собираем уникальные валидные (CHROM, START)
        for _, row in df_filtered[['CHROM', 'START']].iterrows():
            valid_window_starts.add((row['CHROM'], row['START']))

    except Exception as e:
        sys.stderr.write(f"Ошибка при чтении BED файла для фильтрации: {e}\n")
        sys.exit(1)

    return valid_window_starts


def run_window_analysis(df_vcf, scaffold_lengths, exclude_set, bed_filter_set=None):
    """
    Унифицированная логика скользящих окон (1Мб/100кб).
    Если bed_filter_set передан, он используется для отбора генерируемых окон.
    """
    window_size = 1000000
    step = 100000
    kb_divisor = window_size / 1000.0

    # 1. Генерируем полный список ВАЛИДНЫХ окон.
    window_data = []

    for chrom, length in scaffold_lengths.items():
        # Игнорируем исключенные скаффолды и слишком короткие
        if chrom in exclude_set or length < window_size:
            continue

        # Генерируем все потенциальные старты для скользящих окон
        starts = np.arange(0, length - window_size + 1, step)

        for start in starts:
            key = (chrom, start)

            # Применяем фильтр BED, если он есть
            if bed_filter_set is None or key in bed_filter_set:
                 window_data.append((chrom, start, start + window_size))

    if not window_data:
        return []

    # Преобразуем в DataFrame окон
    df_windows = pd.DataFrame(window_data, columns=['CHROM', 'START', 'END'])

    # Инициализируем счетчики
    window_counts = pd.Series(0, index=df_windows.index, dtype=np.int32)
    grouped_vcf = df_vcf.groupby('CHROM')

    # 2. Подсчет SNP в окнах (самая долгая часть, но векторизована внутри цикла)

    for chrom in df_windows['CHROM'].unique():

        try:
            snp_subset = grouped_vcf.get_group(chrom)
        except KeyError:
            continue

        window_subset = df_windows[df_windows['CHROM'] == chrom]

        # Итерация по окнам с векторизованным подсчетом SNP
        for index, row in window_subset.iterrows():
            start = row['START']
            end = row['END']

            # Векторизованный подсчет: сколько SNP попадают в [start, end)
            count = ((snp_subset['POS_0BASED'] >= start) & (snp_subset['POS_0BASED'] < end)).sum()
            window_counts.loc[index] = count

    # 3. Расчет плотности (SNP/kb)
    snp_densities = (window_counts / kb_divisor).tolist()

    return snp_densities


def main():
    parser = argparse.ArgumentParser(description="Анализ плотности SNP с использованием Pandas.")
    parser.add_argument("--vcf", required=True, help="Путь к VCF файлу.")
    parser.add_argument("--scaffolds", required=True,
                        help="Путь к TSV файлу с длинами скаффолдов (scaffold_name\\tlength).")
    parser.add_argument("--bed", required=False, default=None,
                        help="(Опционально) Путь к BED файлу. Если указан, используется для фильтрации скользящих окон по 4-му столбцу.")
    parser.add_argument("--filter_threshold", type=float, default=75000.0, help="Порог для 4-го столбца BED (используется только при указании --bed).")
    parser.add_argument("--exclude", required=False, nargs='+', help="Список скаффолдов/хромосом для исключения.")

    args = parser.parse_args()
    exclude_set = set(args.exclude) if args.exclude else set()

    try:
        # 1. Чтение данных
        scaffold_lengths = parse_scaffold_lengths_tsv(args.scaffolds)
        df_vcf = read_vcf_data(args.vcf, exclude_set)
        total_snps_final = len(df_vcf)

        bed_filter_set = None
        if args.bed:
            # 2. Парсим BED для получения набора (CHROM, START) валидных окон
            bed_filter_set = parse_bed_filter(args.bed, args.filter_threshold)

            if not bed_filter_set:
                sys.stderr.write("Предупреждение: BED файл указан, но ни одно окно не прошло фильтр. Плотности будут 0.0.\n")

        # 3. Запуск унифицированного анализа
        densities = run_window_analysis(df_vcf, scaffold_lengths, exclude_set, bed_filter_set)

        if densities:
            mean_density = statistics.mean(densities)
            median_density = statistics.median(densities)
        else:
            mean_density = 0.0
            median_density = 0.0

        # Формат вывода: Total_SNPs Mean_Density Median_Density
        print(f"{total_snps_final}\t{mean_density:.4f}\t{median_density:.4f}")

    except Exception as e:
        sys.stderr.write(f"Критическая ошибка: {e}\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
