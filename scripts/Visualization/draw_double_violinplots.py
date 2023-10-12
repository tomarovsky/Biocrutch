#!/usr/bin/env python
__author__ = 'tomarovsky'
import argparse
from collections import OrderedDict
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns


def main():
    plt.rcParams.update({'font.size': args.font_size})

    with open(args.input1, "r") as in_fd:
        file_dict = OrderedDict([line.strip().split("\t") for line in in_fd ])

    with open(args.input2, "r") as in_fd:
        file_dict2 = OrderedDict([line.strip().split("\t") for line in in_fd ])

    df_dict = OrderedDict({})
    df_dict2 = OrderedDict({})

    for entry in file_dict:
        df_dict[entry] = pd.read_csv(file_dict[entry], sep="\t", index_col=["CHROM",])
        df_dict[entry]["density"] = df_dict[entry][df_dict[entry].columns[-1]] / args.window_size * args.multiplicator

    for entry in file_dict2:
        df_dict2[entry] = pd.read_csv(file_dict2[entry], sep="\t", index_col=["CHROM",])
        df_dict2[entry]["density"] = df_dict2[entry][df_dict2[entry].columns[-1]] / args.window_size * args.multiplicator

    fig, ax = plt.subplots(figsize=(args.figure_width_per_sample * len(file_dict), args.figure_height), dpi=args.dpi)
    plt.xticks(rotation=args.rotation, ha='right') #, rotation_mode='anchor')

    df_dict_den = OrderedDict({})
    for k,v in df_dict.items():
        df_dict_den[k] = v["density"]

    merge_df = pd.concat(df_dict_den, axis=1)

    df_dict_den2 = OrderedDict({})
    for k,v in df_dict2.items():
        df_dict_den2[k] = v["density"]

    merge_df2 = pd.concat(df_dict_den2, axis=1)

    result_df_list = []
    for col in merge_df.columns:
        for value in merge_df[col]:
            d = {'density': value, 'id': col, 'Reference': f'{args.legend_labels_list[0]}'}
            result_df_list.append(d)
    result_df = pd.DataFrame(result_df_list, columns=('density', 'id', 'Reference'))
    print(result_df)
    result_df_list2 = []
    for col in merge_df2.columns:
        for value in merge_df2[col]:
            d = {'density': value, 'id': col, 'Reference': f'{args.legend_labels_list[1]}'}
            result_df_list2.append(d)
    result_df2 = pd.DataFrame(result_df_list2, columns=('density', 'id', 'Reference'))
    print(result_df2)
    result = pd.concat([result_df, result_df2], axis=0, ignore_index=True, sort=False)
    print(result)

    sns.violinplot(data=result, x="id", y="density", hue="Reference", split = True, scale='width', inner='box', linewidth=1.5, saturation=1, boxprops={'alpha' : 1}, palette=args.colors_list)

    ax.set_xticklabels(list(df_dict.keys()))
    plt.yticks(args.yticklist)
    if args.horizontal_lines:
        for ycoord in args.horizontal_lines:
            plt.axhline(y=ycoord, color="#CC0000", linestyle="--", linewidth=1)
    plt.ylim(ymax=args.ymax, ymin=args.ymin)
    plt.subplots_adjust(left=args.subplots_adjust_left, right=args.subplots_adjust_right,
                        top=args.subplots_adjust_top, bottom=args.subplots_adjust_bottom)
    if args.figure_grid:
        plt.grid(color="gray", linestyle = '--', linewidth = 1.5)
    ax.set_ylabel(args.ylabel)
    ax.set_xlabel('')
    plt.title(args.title)
    plt.legend(loc='upper left', title='Reference')
    plt.tight_layout()
    for ext in args.output_formats:
        plt.savefig("{0}.{1}".format(args.output_prefix, ext))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input1", required=True,
                        help="Input file with two columns containing label in the first one and filename in the second.")
    parser.add_argument("--input2", required=True,
                        help="Input file with two columns containing label in the first one and filename in the second.")
    parser.add_argument("--legend_labels_list", action="store", dest="legend_labels_list", type=lambda s: s.split(","),
                        default=("Species_1", "Species_2"),
                        help="Comma-separated list of species labels of output figure (for input1 and input2, respectively)."
                        "Default: Species_1, Species_2")
    parser.add_argument("--colors_list", action="store", dest="colors_list", type=lambda s: s.split(","),
                        default=('#6094C3', '#E04B4B'),
                        help="Comma-separated list of colors of output figure (for input1 and input2, respectively)."
                        "Default: #6094C3 (blue), #E04B4B (red)")
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
                        default=("svg", "png"),
                        help="Comma-separated list of formats (supported by matlotlib) of "
                             "output figure.Default: svg,png")

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

    args = parser.parse_args()
    print(args.legend_labels_list)
    main()


