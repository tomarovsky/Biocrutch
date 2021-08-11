#!/usr/bin/env python3
__author__ = 'tomarovsky'
from argparse import ArgumentParser

def main():
    result = []

    if len(args.input) == 4:
        with open(args.input[0], "r") as f:
            for i in f:
                if len(i) > 2:
                    axis_x, axis_y = i.split()[0].replace(".", ","), i.split()[1].replace(".", ",")
                    tmp = [axis_x, axis_y, "", "", ""]
                    result.append(tmp)
                else:
                    result.append("")
        with open(args.input[1], "r") as f:
            for i in f:
                if len(i) > 2:
                    axis_x, axis_y = i.split()[0].replace(".", ","), i.split()[1].replace(".", ",")
                    tmp = [axis_x, "", axis_y, "", ""]
                    result.append(tmp)
                else:
                    result.append("")
        with open(args.input[2], "r") as f:
            for i in f:
                if len(i) > 2:
                    axis_x, axis_y = i.split()[0].replace(".", ","), i.split()[1].replace(".", ",")
                    tmp = [axis_x, "", "", axis_y, ""]
                    result.append(tmp)
                else:
                    result.append("")
        with open(args.input[3], "r") as f:
            for i in f:
                if len(i) > 2:
                    axis_x, axis_y = i.split()[0].replace(".", ","), i.split()[1].replace(".", ",")
                    tmp = [axis_x, "", "", "", axis_y]
                    result.append(tmp)
                else:
                    result.append("")
    elif len(args.input) == 2:
        with open(args.input[0], "r") as f:
            for i in f:
                if len(i) > 2:
                    axis_x, axis_y = i.split()[0].replace(".", ","), i.split()[1].replace(".", ",")
                    tmp = [axis_x, axis_y, "", "", ""]
                    result.append(tmp)
                else:
                    result.append("")
        with open(args.input[1], "r") as f:
            for i in f:
                if len(i) > 2:
                    axis_x, axis_y = i.split()[0].replace(".", ","), i.split()[1].replace(".", ",")
                    tmp = [axis_x, "", "", axis_y, ""]
                    result.append(tmp)
                else:
                    result.append("")

    with open(args.output, "a") as res:
        for i in result:
            line = "\t".join(i)
            # print(line)
            res.write(line)
            res.write("\n")




if __name__ == "__main__":
    parser = ArgumentParser(description="script for combining PSMC data into one table for visualization")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str,
                                nargs='+', help="PSMC files (first rounds(1 and 2), then diploids(1 and 2)!)")
    group_required.add_argument('-o', '--output', type=str, help="outfile name", default="PSMC.results.tab")
    args = parser.parse_args()
    main()
