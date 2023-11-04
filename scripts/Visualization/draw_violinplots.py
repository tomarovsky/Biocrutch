#!/usr/bin/env python
__author__ = 'tomarovsky'
from collections import OrderedDict
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import argparse


def main():
    if args.pretty:
        sns.set_theme(style="ticks")
        sns.set_context("talk")
    else:
        plt.rcParams.update({'font.size': args.font_size})

    with open(args.input, "r") as in_fd:
        file_dict = OrderedDict([line.strip().split("\t") for line in in_fd ])

    df_dict = OrderedDict({})

    for entry in file_dict:
        df_dict[entry] = pd.read_csv(file_dict[entry], sep="\t", index_col=["CHROM",])
        print(df_dict[entry])
        df_dict[entry]["density"] = df_dict[entry][df_dict[entry].columns[1]] / args.window_size * args.multiplicator

    fig, ax = plt.subplots(figsize=(args.figure_width_per_sample * len(file_dict), args.figure_height), dpi=args.dpi)
    plt.xticks(rotation=args.rotation)

    sns.violinplot(data=[df_dict[entry]["density"] for entry in df_dict], palette='turbo', scale='width',
                   inner=None, linewidth=0, saturation=0.4) # or inner="box"
    sns.boxplot(data=[df_dict[entry]["density"] for entry in df_dict], palette='turbo', width=0.1,
            linewidth=0.4, flierprops={'marker':'o', 'markersize':1}, boxprops={'zorder': 2}, ax=ax)

    ax.set_xticklabels(list(df_dict.keys()), ha='right', rotation_mode='anchor')
    plt.yticks(args.yticklist)
    if args.horizontal_lines:
        for ycoord in args.horizontal_lines:
            plt.axhline(y=ycoord, color="red", linestyle="--", linewidth=1)
    plt.ylim(ymax=args.ymax, ymin=args.ymin)
    plt.subplots_adjust(left=args.subplots_adjust_left, right=args.subplots_adjust_right,
                        top=args.subplots_adjust_top, bottom=args.subplots_adjust_bottom)
    if args.figure_grid:
        plt.grid(color="gray", linestyle = '--', linewidth = 0.5)
    ax.set_ylabel(args.ylabel)
    ax.set_xlabel("Sample name")
    plt.title(args.title)
    plt.tight_layout()
    for ext in args.output_formats:
        plt.savefig("{0}.{1}".format(args.output_prefix, ext))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True,
                        help="Input file with two columns containing label in the first one and filename in the second."
                             "Boxplots will be drawn in the same order as labels")
    parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                        help="Prefix of output files")
    parser.add_argument("-d", "--dpi", action="store", dest="dpi", type=int, default=300,
                        help="Dpi of figure")

    parser.add_argument("--figure_height", action="store", dest="figure_height", type=float, default=6.0,
                        help="Height of figure in inches. Default: 6")
    parser.add_argument("--figure_width_per_sample", action="store", dest="figure_width_per_sample", type=float, default=1,
                        help="Per sample width of figure in inches. Default: 1")
    parser.add_argument("--figure_grid", action="store_true", default=False,
                        help="Add grid lines to the figure. Default: False")
    parser.add_argument("--font-size", action="store", dest="font_size", type=float, default=16,
                        help="Font size. Default: 16")
    parser.add_argument("-e", "--output_formats", action="store", dest="output_formats", type=lambda s: s.split(","),
                        default=("png", "svg"), help="Comma-separated list of formats (supported by matlotlib) of "
                             "output figure.Default: png,svg")

    parser.add_argument("-l", "--title", action="store", dest="title", default="Variant density",
                        help="Suptitle of figure. Default: 'Variant density'")
    parser.add_argument("--ylabel", action="store", dest="ylabel", default="Heterozygous SNPs/kbp",
                                help="y-axis label. Default: 'Heterozygous SNPs/kbp'")
    parser.add_argument("-w", "--window_size", action="store", dest="window_size", required=True, type=float,
                        help="Size of the windows use for counts.")
    parser.add_argument("-m", "--multiplicator", action="store", dest="multiplicator", default=1000, type=float,
                        help="Multiplicator for variant counts. "
                             "Default: 1000, i.e variant counts will be scaled to per 1 kbp ")
    parser.add_argument("--ymin", action="store", dest="ymin", type=float, default=-0.1,
                        help="Minimum limit for Y axis . Default: -0.1")
    parser.add_argument("--ymax",  action="store", dest="ymax", type=float, default=None,
                        help="Maximum limit for Y axis. Default: not set")
    parser.add_argument("--yticklist",  action="store", dest="yticklist", type=lambda s: list(map(float, s.split(","))),
                        default=[0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 2, 3, 4, 5],
                        help="Comma-separated tick list for Y axis. "
                             "Default: 0.05,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,2,3,4,5")
    parser.add_argument("--rotation",  action="store", dest="rotation", type=float, default=90,
                        help="Rotation angle for X labels. Default: 90")
    parser.add_argument("--horizontal_lines",  action="store", dest="horizontal_lines",
                        type=lambda s: list(map(float, s.split(","))),
                        help="Comma-separated list of y-coordinates to draw horizontal lines. "
                             "Default: not set")

    parser.add_argument("--subplots_adjust_left", action="store", dest="subplots_adjust_left", type=float,
                        help="Adjust left border of subplots on the figure. Default: matplotlib defaults")
    parser.add_argument("--subplots_adjust_top", action="store", dest="subplots_adjust_top", type=float,
                        help="Adjust top border of subplots on the figure. Default: matplotlib defaults")
    parser.add_argument("--subplots_adjust_right", action="store", dest="subplots_adjust_right", type=float,
                        help="Adjust right border of subplots on the figure. Default: matplotlib defaults")
    parser.add_argument("--subplots_adjust_bottom", action="store", dest="subplots_adjust_bottom", type=float,
                        help="Adjust bottom border of subplots on the figure. Default: matplotlib defaults")
    parser.add_argument("--only_count", action="store_true", dest="only_count", default=False,
                        help="Only count variants, do not draw them. Default: False")
    parser.add_argument("--pretty", action="store_true", dest="pretty", default=False,
                        help="Pretty visualization using theme 'ticks' and context 'talk' from Seaborn. Default: False")

    args = parser.parse_args()
    main()


