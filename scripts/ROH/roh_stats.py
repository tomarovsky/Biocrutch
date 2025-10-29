#!/usr/bin/env python3
__author__ = 'tomarovsky'
import argparse

import numpy as np
import pandas as pd


def analyze_roh(input_file: str, genome_length: int, exclude_chr: str = None):
    df = pd.read_csv(
        input_file,
        sep='\s+',  # Разделитель: пробелы (для вашего примера данных)
        header=None,
        names=['Scaffold', 'Start', 'End', 'Length']
    )

    if exclude_chr:
        exclude_list = [c.strip() for c in exclude_chr.split(',') if c.strip()]
        if exclude_list:
            df = df[~df['Scaffold'].isin(exclude_list)].copy()

    Number_of_ROHs = len(df)
    Total_length = df['Length'].sum()
    Percent_of_ROHs = (Total_length / genome_length) * 100
    # Classification of ROHs by length
    df['Class'] = np.select(
        [
            df['Length'] < 1_000_000,
            (df['Length'] >= 1_000_000) & (df['Length'] < 10_000_000),
            df['Length'] >= 10_000_000
        ],
        [
            'S_ROHS',
            'L_ROHS',
            'UL_ROHS'
        ],
        default='N/A'
    )

    class_stats = df.groupby('Class')['Length'].agg(['count', 'sum']).reindex(['S_ROHS', 'L_ROHS', 'UL_ROHS']).fillna(0)

    # 4. Number_of_S_ROHS и 5. Percent_of_S_ROHs
    Number_of_S_ROH = int(class_stats.loc['S_ROHS', 'count'])
    Length_of_S_ROH = class_stats.loc['S_ROHS', 'sum']
    Percent_of_S_ROH = (Length_of_S_ROH / genome_length) * 100

    # 6. Number_of_L_ROHS и 7. Percent_of_L_ROHs
    Number_of_L_ROH = int(class_stats.loc['L_ROHS', 'count'])
    Length_of_L_ROH = class_stats.loc['L_ROHS', 'sum']
    Percent_of_L_ROH = (Length_of_L_ROH / genome_length) * 100

    # 8. Number_of_UL_ROHS и 9. Percent_of_UL_ROHS
    Number_of_UL_ROH = int(class_stats.loc['UL_ROHS', 'count'])
    Length_of_UL_ROH = class_stats.loc['UL_ROHS', 'sum']
    Percent_of_UL_ROH = (Length_of_UL_ROH / genome_length) * 100

    # Формирование списка результатов
    result_metrics = [
        Number_of_ROHs, Total_length / 1000000, Percent_of_ROHs,
        Number_of_S_ROH, Percent_of_S_ROH,
        Number_of_L_ROH, Percent_of_L_ROH,
        Number_of_UL_ROH, Percent_of_UL_ROH
    ]

    # --- Форматирование и вывод ---

    # Заголовок
    header = [
        "Number_of_ROHs", "Total_length", "Percent_of_ROHs",
        "Number_of_S_ROH", "Percent_of_S_ROH",
        "Number_of_L_ROH", "Percent_of_L_ROH",
        "Number_of_UL_ROH", "Percent_of_UL_ROH"
    ]

    # Форматирование значений
    output_values = [
        f"{result_metrics[0]:d}",        # Number_of_ROHs (целое)
        f"{result_metrics[1]:.2f}",      # Total_length (2 знаков после запятой)
        f"{result_metrics[2]:.2f}",      # Percent_of_ROHs (2 знаков после запятой)
        f"{result_metrics[3]:d}",        # Number_of_S_ROH (целое)
        f"{result_metrics[4]:.2f}",      # Percent_of_S_ROH (2 знаков после запятой)
        f"{result_metrics[5]:d}",        # Number_of_L_ROH (целое)
        f"{result_metrics[6]:.2f}",      # Percent_of_L_ROH (2 знаков после запятой)
        f"{result_metrics[7]:d}",        # Number_of_UL_ROH (целое)
        f"{result_metrics[8]:.2f}"       # Percent_of_UL_ROH (2 знаков после запятой)
    ]

    # Вывод заголовка и данных
    # print('\t'.join(header))
    print('\t'.join(output_values))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="ROH statistics",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '-i', '--input_file',
        type=str,
        help="path to .roh file from Biocrutch"
    )
    parser.add_argument(
        '-g', '--genome_length',
        type=int,
        help="total genome length"
    )
    parser.add_argument(
        '-e', '--exclude_chr',
        type=str,
        default=None,
        help="exclude scaffolds, separated by commas (e.g.: chrX)"
    )

    args = parser.parse_args()

    analyze_roh(args.input_file, args.genome_length, args.exclude_chr)
