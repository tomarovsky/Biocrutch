#!/usr/bin/env python3
__author__ = 'tomarovsky'
from Biocrutch.Routines.routine_functions import metaopen, metaoutput
import argparse
from sys import stdin

def main():
    if args.input == stdin:
        outprefix = "repeatmasker.out"
    else:
        if args.output:
            outprefix = args.output
        else:
            outprefix = ('.').join(args.input.split('.')[:-1])
    outprefix = metaoutput(outprefix, ".gff.gz")
    outfile = metaopen(outprefix, "wt")
    print ("OUT to GFF converter...")
    with metaopen(args.input, "r", buffering=args.buffering) as data:
        count = 0
        for l in data:
            if count < 4:
                count += 1
                continue
            line = l.strip().split()
            strand = "+" if line[8] == "+" else "-"
            class_family_info = line[10].split("/")
            if len(class_family_info) == 1:
                class_family_info.append(".")
            additional_info = "Class=%s;Family=%s;Matching_repeat=%s;SW_score=%s;Perc_div=%s;Perc_del=%s;Pers_ins=%s"\
                               % (class_family_info[0], class_family_info[1],
                                  line[9], line[0], line[1], line[2], line[3])
            gff_line = "%s\tRepeatMasker\trepeat\t%s\t%s\t.\t%s\t.\t%s\n" % (line[4], line[5], line[6], strand, additional_info)
            outfile.write(gff_line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="repeatmasker.out to GFF converter")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str, default=stdin,
                                help="input file or stdin")
    group_additional = parser.add_argument_group('Additional options')
    group_additional.add_argument('-b', '--buffering', type=int, 
                                  default=None, help="Text buffering. Default = None")
    group_additional.add_argument('-o', '--output', type=str,
                                  help="output file prefix")
    args = parser.parse_args()
    main()