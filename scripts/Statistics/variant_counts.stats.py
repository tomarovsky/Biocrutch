#!/usr/bin/env python3
__author__ = 'tomarovsky'

import argparse
import sys
import statistics
import gzip
import pandas as pd
from io import StringIO
import numpy as np


def open_file(filename):
    """Возвращает файловый дескриптор или gzip.open."""
    if filename.endswith(".gz"):
        return gzip.open(filename, 'rt', encoding='utf-8')
    else:
        return open(filename, 'r', encoding='utf-8')

# --- Функции парсинга ---

def parse_scaffold_lengths_tsv(tsv_file):
    """Парсит TSV файл с длинами скаффолдов с помощью Pandas."""
    try:
        # Читаем TSV, используя первый столбец как индекс (имя скаффолда)
        df = pd.read_csv(tsv_file, sep='\s+', header=None, names=['CHROM', 'LENGTH'], dtype={'CHROM': str, 'LENGTH': np.int64})

        if df.empty:
            raise ValueError(f"Файл {tsv_file} не содержит валидных данных.")

        # Преобразуем DataFrame в словарь {scaffold: length}
        scaffold_lengths = df.set_index('CHROM')['LENGTH'].to_dict()
        return scaffold_lengths

    except Exception as e:
        sys.stderr.write(f"Ошибка при чтении TSV с длинами: {e}\n")
        sys.exit(1)


def read_vcf_data(vcf_file, exclude_set):
    """
    Читает VCF файл (только данные, без заголовка) в Pandas DataFrame.
    Сразу фильтрует по exclude_set.
    """
    header_lines = 0

    # 1. Находим, где заканчивается заголовок VCF
    with open_file(vcf_file) as f:
        for line in f:
            if line.startswith("#"):
                header_lines += 1
            else:
                break

    # 2. Читаем данные VCF
    try:
        df_vcf = pd.read_csv(
            vcf_file,
            sep='\t',
            skiprows=header_lines - 1, # Пропускаем строки заголовка, кроме последней (#CHROM...)
            header=0,
            usecols=[0, 1], # Интересуют только CHROM и POS
            names=['CHROM', 'POS'],
            dtype={'CHROM': str, 'POS': np.int32}
        )
        # Отсчет позиции делаем 0-based
        df_vcf['POS_0BASED'] = df_vcf['POS'] - 1

        # 3. Фильтрация по exclude_set
        total_snps_before_exclude = len(df_vcf)
        if exclude_set:
             df_vcf = df_vcf[~df_vcf['CHROM'].isin(exclude_set)]

        return df_vcf, total_snps_before_exclude

    except Exception as e:
        sys.stderr.write(f"Ошибка при чтении VCF данных: {e}\n")
        sys.exit(1)

# --- Логика анализа ---

def run_with_bed_pandas(df_vcf, bed_file, threshold, exclude_set):
    """
    Логика VCF + BED с использованием Pandas.
    Это самая быстрая часть.
    """

    # 1. Читаем и фильтруем BED
    try:
        # Читаем только 4 столбца: CHROM, START, END, VALUE
        df_bed = pd.read_csv(
            bed_file,
            sep='\t',
            comment='#',
            header=None,
            usecols=[0, 1, 2, 3],
            names=['CHROM', 'START', 'END', 'VALUE'],
            dtype={'CHROM': str, 'START': np.int32, 'END': np.int32, 'VALUE': np.float64}
        )
    except Exception as e:
        sys.stderr.write(f"Ошибка при чтении BED файла: {e}\n")
        sys.exit(1)

    # Применяем фильтры
    df_bed = df_bed[~df_bed['CHROM'].isin(exclude_set)]
    df_bed = df_bed[df_bed['VALUE'] <= threshold]

    if df_bed.empty:
         return len(df_vcf), []

    # 2. Определение окна и kb_divisor
    # Берем размер первого валидного окна как kb_divisor
    window_size = (df_bed['END'].iloc[0] - df_bed['START'].iloc[0])
    kb_divisor = window_size / 1000.0

    # 3. Присвоение SNP окнам (Spatial Join с Pandas/NumPy)

    # Создаем индекс для быстрого поиска
    df_bed = df_bed.set_index(['CHROM', 'START', 'END'])

    # Инициализируем счетчики
    window_counts = pd.Series(0, index=df_bed.index, dtype=np.int32)

    # Перебираем скаффолды, присутствующие в валидных окнах BED
    for chrom in df_bed.index.get_level_values('CHROM').unique():
        # Фильтруем SNP по текущему скаффолду
        snp_subset = df_vcf[df_vcf['CHROM'] == chrom]

        # Фильтруем окна по текущему скаффолду
        bed_subset = df_bed.xs(chrom, level='CHROM')

        # Основной цикл по окнам (трудно векторизовать без специализированных либ)
        for idx in bed_subset.index:
            start, end = idx

            # Векторизованный подсчет SNP, попадающих в текущее окно [start, end)
            count = ((snp_subset['POS_0BASED'] >= start) & (snp_subset['POS_0BASED'] < end)).sum()
            window_counts.loc[(chrom, start, end)] = count

    # 4. Расчет плотности
    snp_densities = (window_counts / kb_divisor).tolist()

    return len(df_vcf), snp_densities


def run_vcf_only_pandas(df_vcf, scaffold_lengths, exclude_set):
    """
    Логика VCF-only (скользящие окна) с использованием Pandas.
    """
    window_size = 1000000
    step = 100000
    kb_divisor = window_size / 1000.0

    # 1. Генерируем список ВАЛИДНЫХ окон.
    window_data = []

    for chrom, length in scaffold_lengths.items():
        if chrom in exclude_set or length < window_size:
            continue

        # Создаем окна, которые ПОЛНОСТЬЮ помещаются (векторизованно)
        starts = np.arange(0, length - window_size + 1, step)
        ends = starts + window_size

        # Добавляем данные в список
        for start, end in zip(starts, ends):
            window_data.append((chrom, start, end))

    if not window_data:
        return len(df_vcf), []

    # Преобразуем в DataFrame окон
    df_windows = pd.DataFrame(window_data, columns=['CHROM', 'START', 'END'])
    df_windows['COUNT'] = 0

    # 2. Подсчет SNP в окнах

    # Инициализируем счетчики (будем использовать их как Series)
    window_counts = pd.Series(0, index=df_windows.index, dtype=np.int32)

    # Группируем SNP по скаффолдам
    grouped_vcf = df_vcf.groupby('CHROM')

    # Итерация по скаффолдам, где были сгенерированы окна
    for chrom in df_windows['CHROM'].unique():

        try:
            snp_subset = grouped_vcf.get_group(chrom)
        except KeyError:
            # Если скаффолд есть в --scaffolds, но нет SNP в VCF, пропускаем
            continue

        # Фильтруем окна для текущего скаффолда
        window_subset = df_windows[df_windows['CHROM'] == chrom]

        # --- Векторизованная логика подсчета перекрытий ---
        # В этом режиме одна SNP может попасть в 10 окон (1Мб/100кб),
        # поэтому нам нужно найти все перекрытия.

        for index, row in window_subset.iterrows():
            start = row['START']
            end = row['END']

            # Векторизованный подсчет SNP, попадающих в текущее окно [start, end)
            count = ((snp_subset['POS_0BASED'] >= start) & (snp_subset['POS_0BASED'] < end)).sum()
            window_counts.loc[index] = count

    # 3. Расчет плотности
    snp_densities = (window_counts / kb_divisor).tolist()

    return len(df_vcf), snp_densities


def main():
    parser = argparse.ArgumentParser(description="Анализ плотности SNP с использованием Pandas.")
    parser.add_argument("--vcf", required=True, help="Путь к VCF файлу.")
    parser.add_argument("--scaffolds", required=True,
                        help="Путь к TSV файлу с длинами скаффолдов (scaffold_name\\tlength).")
    parser.add_argument("--bed", required=False, default=None,
                        help="(Опционально) Путь к BED файлу с окнами. Если не указан, используется режим скользящих окон (1Мб/100кб).")
    parser.add_argument("--filter_threshold", type=float, default=75000.0, help="Порог для 4-го столбца BED (игнорируется в VCF-only).")
    parser.add_argument("--exclude", required=False, nargs='+', help="Список скаффолдов/хромосом для исключения.")

    args = parser.parse_args()
    exclude_set = set(args.exclude) if args.exclude else set()

    total_snps = 0
    densities = []

    try:
        # 1. Чтение длин скаффолдов
        scaffold_lengths = parse_scaffold_lengths_tsv(args.scaffolds)

        # 2. Чтение и предварительная фильтрация VCF
        df_vcf, total_snps_with_exclude = read_vcf_data(args.vcf, exclude_set)

        # total_snps_final - это количество SNP, которые не попали в --exclude
        total_snps_final = len(df_vcf)

        if total_snps_final == 0:
            sys.stderr.write("Предупреждение: После фильтрации --exclude не осталось SNP для анализа.\n")

        if args.bed:
            total_snps, densities = run_with_bed_pandas(
                df_vcf, args.bed, args.filter_threshold, exclude_set
            )
        else:
            total_snps, densities = run_vcf_only_pandas(df_vcf, scaffold_lengths, exclude_set)

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
