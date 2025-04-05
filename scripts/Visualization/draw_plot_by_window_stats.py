#!/usr/bin/env python3
__author__ = 'tomarovsky'
import numpy as np
import matplotlib.pyplot as plt

plt.ioff()
from argparse import ArgumentParser


def draw_plot_by_window_stats(
    input_file,
    output_prefix,
    metric,
    separator='\t',
    min_x=None,
    max_x=None,
    min_y=None,
    max_y=None,
    extensions=['png', 'svg'],
    xlabel=None,
    ylabel=None,
    title=None,
    width=6,
    height=6,
    markersize=2,
    ylogbase=10,
    type='plot',
    grid=False,
    close_plot=True,
):
    if metric.lower() == 'median':
        data = np.loadtxt(input_file, comments='#', usecols=(1, 2), delimiter=separator)
    elif metric.lower() == 'average':
        data = np.loadtxt(input_file, comments='#', usecols=(1, 3), delimiter=separator)
    elif metric.lower() == 'max':
        data = np.loadtxt(input_file, comments='#', usecols=(1, 4), delimiter=separator)
    elif metric.lower() == 'min':
        data = np.loadtxt(input_file, comments='#', usecols=(1, 5), delimiter=separator)

    plt.figure(1, figsize=(width, height), dpi=300)
    plt.subplot(1, 1, 1)
    if type == 'plot':
        plt.plot(data[:, 0], data[:, 1], markersize=markersize)
    elif type == 'scatter':
        plt.scatter(data[:, 0], data[:, 1], s=markersize)
    plt.xlim(xmin=min_x, xmax=max_x)
    plt.ylim(ymin=min_y, ymax=max_y)
    if xlabel:
        plt.xlabel(xlabel)
    if ylabel:
        plt.ylabel(ylabel)
    if title:
        plt.title(title)
    if grid or grid == 'True':
        plt.grid()
    for ext in extensions:
        plt.savefig(f'{output_prefix}.{type}.{ext}')

    plt.yscale('log')

    for ext in extensions:
        plt.savefig(f'{output_prefix}.{type}.{ylogbase}.{ext}')
    if close_plot:
        plt.close()


def main():
    draw_plot_by_window_stats(
        args.input_file,
        args.output_prefix,
        args.metric,
        separator=args.separator,
        extensions=args.extensions,
        min_x=args.min_x,
        max_x=args.max_x,
        min_y=args.min_y,
        max_y=args.max_y,
        xlabel=args.xlabel,
        ylabel=args.ylabel,
        title=args.title,
        width=args.width,
        height=args.height,
        markersize=args.markersize,
        ylogbase=args.ylogbase,
        type=args.type,
        grid=args.grid,
        close_plot=args.close_plot,
    )


if __name__ == '__main__':
    parser = ArgumentParser(description='coverage visualization from stats.csv file (scripts/Coverage/coverage_statistics.py)')
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input-file', type=str, help='mosdepth.bed.gz file')
    group_required.add_argument('-o', '--output-prefix', type=str, help='output prefix')
    group_additional = parser.add_argument_group('Additional options')
    group_additional.add_argument('-m', '--metric', type=str, help="metric: 'median', 'average', 'max' or 'min'", default='median')
    group_additional.add_argument('-s', '--separator', type=str, help='separator', default='\t')
    group_additional.add_argument(
        '-e', '--extensions', type=lambda s: [str(ext) for ext in s.split(',')], default=['png', 'svg'], help='output files extensions'
    )
    group_additional.add_argument('--min_x', help='min_x value', type=int, default=None)
    group_additional.add_argument('--max_x', help='max_x value', type=int, default=None)
    group_additional.add_argument('--min_y', help='min_y value', type=int, default=None)
    group_additional.add_argument('--max_y', help='max_y value', type=int, default=None)
    group_additional.add_argument('--xlabel', help='xlabel', default=None)
    group_additional.add_argument('--ylabel', help='ylabel', default=None)
    group_additional.add_argument('--title', type=str, help='title', default=None)
    group_additional.add_argument('--width', type=int, help='xlabel', default=6)
    group_additional.add_argument('--height', type=int, help='ylabel', default=6)
    group_additional.add_argument('--markersize', type=int, help='markersize', default=2)
    group_additional.add_argument('--ylogbase', type=int, help='ylogbase', default=10)
    group_additional.add_argument('--type', type=str, help='type', default='plot')
    group_additional.add_argument('--grid', type=bool, help='grid', default=False)
    group_additional.add_argument('--close_plot', type=bool, help='close plot', default=True)
    args = parser.parse_args()
    main()
