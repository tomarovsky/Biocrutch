#!/usr/bin/env python3
__author__ = 'tomarovsky'
from argparse import ArgumentParser

def main():
    result = []

    number_of_files = len(args.input)

    c = 1
    flag = False
    for i in range(0, number_of_files):
        with open(args.input[i], "r") as f:
            for line in f:
                if len(line) > 2:
                    axis_x, axis_y = line.split()[0].replace(".", ","), line.split()[1].replace(".", ",")
                    tmp = [axis_x] + [""] * 10
                   
                    if i < number_of_files/2:
                        tmp[i + 1] = axis_y
                        print(tmp)
                    else:
                        tmp[int(5 + c)] = axis_y
                        flag = True
                        print(tmp)
                    result.append(tmp)
                else:
                    result.append("")
        if flag == True:
            c += 1

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
