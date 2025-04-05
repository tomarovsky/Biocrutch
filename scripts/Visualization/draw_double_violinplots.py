#!/usr/bin/env python
__author__ = 'tomarovsky'
import argparse
from collections import OrderedDict
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def load_and_prepare_data(file_dict, label, args):
    data_list = []
    for entry in file_dict:
        df = pd.read_csv(file_dict[entry], sep='\t')
        if args.no_x:
            df = df[~df['CHROM'].str.contains('chrX')]
        df.set_index('CHROM', inplace=True)
        density = df.iloc[:, -1] / args.window_size * args.multiplicator
        data_list.extend([{'density': d, 'id': entry, 'Reference': label} for d in density])
    return data_list


def main():
    plt.rcParams.update({'font.size': args.font_size})

    with open(args.input1, 'r') as in_fd:
        file_dict1 = OrderedDict([line.strip().split('\t') for line in in_fd])

    with open(args.input2, 'r') as in_fd:
        file_dict2 = OrderedDict([line.strip().split('\t') for line in in_fd])

    data1 = load_and_prepare_data(file_dict1, args.legend_labels_list[0], args)
    data2 = load_and_prepare_data(file_dict2, args.legend_labels_list[1], args)

    result = pd.DataFrame(data1 + data2)
    if args.only_count:
        result.to_csv('counts.csv', index=False)

    fig, ax = plt.subplots(figsize=(args.figure_width_per_sample * len(file_dict1), args.figure_height), dpi=args.dpi)
    plt.xticks(rotation=args.rotation, ha='right')

    sns.violinplot(
        data=result,
        x='id',
        y='density',
        hue='Reference',
        split=True,
        density_norm='width',
        inner='box',
        inner_kws=dict(box_width=1.5, whis_width=0, marker='.', markersize=2.2, markeredgewidth=0.6),
        linewidth=0.5,
        saturation=1,
        palette=args.colors_list,
    )

    ax.set_xticks(range(len(file_dict1.keys())))
    ax.set_xticklabels(file_dict1.keys())
    plt.yticks(args.yticklist)
    if args.horizontal_lines:
        for ycoord in args.horizontal_lines:
            plt.axhline(y=ycoord, color='#CC0000', linestyle='--', linewidth=1)
    plt.ylim(ymax=args.ymax, ymin=args.ymin)
    plt.subplots_adjust(
        left=args.subplots_adjust_left, right=args.subplots_adjust_right, top=args.subplots_adjust_top, bottom=args.subplots_adjust_bottom
    )
    if args.figure_grid:
        plt.grid(color='gray', linestyle='--', linewidth=1.5)
    ax.set_ylabel(args.ylabel)
    ax.set_xlabel('')
    plt.title(args.title)
    plt.legend(loc='upper left', title='Reference', ncol=2)
    plt.tight_layout()
    for ext in args.output_formats:
        plt.savefig('{0}.{1}'.format(args.output_prefix, ext))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--input1', required=True, help='Input file with two columns containing label in the first one and filename in the second.'
    )
    parser.add_argument('--input2', required=True, help='Input file with two columns containing label in the first one and filename in the second.')
    parser.add_argument(
        '--legend_labels_list',
        action='store',
        dest='legend_labels_list',
        type=lambda s: s.split(','),
        default=('Species_1', 'Species_2'),
        help='Comma-separated list of species labels of output figure (for input1 and input2, respectively).Default: Species_1, Species_2',
    )
    parser.add_argument(
        '--colors_list',
        action='store',
        dest='colors_list',
        type=lambda s: s.split(','),
        default=('#6094C3', '#E04B4B'),
        help='Comma-separated list of colors of output figure (for input1 and input2, respectively).Default: #6094C3 (blue), #E04B4B (red)',
    )
    parser.add_argument('-o', '--output_prefix', action='store', dest='output_prefix', required=True, help='Prefix of output files')
    parser.add_argument('-d', '--dpi', action='store', dest='dpi', type=int, default=300, help='Dpi of figure')

    parser.add_argument(
        '--figure_height', action='store', dest='figure_height', type=float, default=6.0, help='Height of figure in inches. Default: 6'
    )
    parser.add_argument(
        '--figure_width_per_sample',
        action='store',
        dest='figure_width_per_sample',
        type=float,
        default=0.5,
        help='Per sample width of figure in inches. Default: 0.5',
    )
    parser.add_argument('--figure_grid', action='store_true', default=False, help='Add grid lines to the figure. Default: False')
    parser.add_argument('--font-size', action='store', dest='font_size', type=float, default=16, help='Font size. Default: 16')
    parser.add_argument(
        '-e',
        '--output_formats',
        action='store',
        dest='output_formats',
        type=lambda s: s.split(','),
        default=('svg', 'png'),
        help='Comma-separated list of formats (supported by matlotlib) of output figure.Default: svg,png',
    )

    parser.add_argument(
        '-l', '--title', action='store', dest='title', default='Variant density', help="Suptitle of figure. Default: 'Variant density'"
    )
    parser.add_argument(
        '--ylabel', action='store', dest='ylabel', default='Heterozygous SNPs/kbp', help="y-axis label. Default: 'Heterozygous SNPs/kbp'"
    )
    parser.add_argument(
        '-w', '--window_size', action='store', dest='window_size', required=True, type=float, help='Size of the windows use for counts.'
    )
    parser.add_argument(
        '-m',
        '--multiplicator',
        action='store',
        dest='multiplicator',
        default=1000,
        type=float,
        help='Multiplicator for variant counts. Default: 1000, i.e variant counts will be scaled to per 1 kbp ',
    )
    parser.add_argument('--ymin', action='store', dest='ymin', type=float, default=-0.1, help='Minimum limit for Y axis . Default: -0.1')
    parser.add_argument('--ymax', action='store', dest='ymax', type=float, default=None, help='Maximum limit for Y axis. Default: not set')
    parser.add_argument(
        '--yticklist',
        action='store',
        dest='yticklist',
        type=lambda s: list(map(float, s.split(','))),
        default=[0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 2, 3, 4, 5],
        help='Comma-separated tick list for Y axis. Default: 0.05,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,2,3,4,5',
    )
    parser.add_argument('--rotation', action='store', dest='rotation', type=float, default=90, help='Rotation angle for X labels. Default: 90')

    parser.add_argument(
        '--horizontal_lines',
        action='store',
        dest='horizontal_lines',
        type=lambda s: list(map(float, s.split(','))),
        help='Comma-separated list of y-coordinates to draw horizontal lines. Default: not set',
    )

    parser.add_argument(
        '--subplots_adjust_left',
        action='store',
        dest='subplots_adjust_left',
        type=float,
        help='Adjust left border of subplots on the figure. Default: matplotlib defaults',
    )
    parser.add_argument(
        '--subplots_adjust_top',
        action='store',
        dest='subplots_adjust_top',
        type=float,
        help='Adjust top border of subplots on the figure. Default: matplotlib defaults',
    )
    parser.add_argument(
        '--subplots_adjust_right',
        action='store',
        dest='subplots_adjust_right',
        type=float,
        help='Adjust right border of subplots on the figure. Default: matplotlib defaults',
    )
    parser.add_argument(
        '--subplots_adjust_bottom',
        action='store',
        dest='subplots_adjust_bottom',
        type=float,
        help='Adjust bottom border of subplots on the figure. Default: matplotlib defaults',
    )
    parser.add_argument('--no_x', action='store_true', dest='no_x', default=False, help='Do not use counts from X chromosome. Default: False')
    parser.add_argument(
        '--only_count',
        action='store_true',
        dest='only_count',
        default=False,
        help='Only count and save variants to CSV, do not draw them. Default: False',
    )

    args = parser.parse_args()

    main()
