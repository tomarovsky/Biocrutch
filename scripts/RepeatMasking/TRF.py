#!/usr/bin/env python3
__author__ = 'tomarovsky'
from Biocrutch.Routines.routine_functions import metaopen, metaoutput
import argparse
from sys import stdin

def main():
    if args.input == stdin:
        outprefix = "trf.dat"
    else:
        if args.output:
            outprefix = args.output
        else:
            outprefix = ('.').join(args.input.split('.')[:-1])
    outfile = metaopen (outprefix, "wt")
    outfile = metaoutput(outfile, ".gff.gz")
    print ("DAT to GFF converter...")
    with metaopen(args.input, "r", buffering=args.buffering) as data:
        count = 0
        for l in data:
            line = l.strip().split(" ")
            if l.startswith('Sequence:'): # information about authors and parameters is not processed
                seq_name = line[1]
            elif l[0].isdigit():
                [start, stop, period, copies, consensus_size, 
                 perc_match, perc_indels, align_score, 
                 perc_A, perc_C, perc_G, perc_T, 
                 entropy, cons_seq, repeat_seq] = line
                gff_line = [seq_name, 'TRF', 'repeat',
                            start, stop, '.', '.', '.', 
                            'ID='+ str(count) + 
                            ';period=' + period +
                            ';copies=' + copies +
                            ';consensus_size=' + consensus_size +
                            ';perc_match=' + perc_match +
                            ';perc_indels=' + perc_indels +
                            ';align_score=' + align_score +
                            'perc_A' + perc_A + 
                            'perc_C' + perc_C + 
                            'perc_G' + perc_G + 
                            'perc_T' + perc_T + 
                            ';entropy=' + entropy +
                            ';cons_seq=' + cons_seq +
                            ';repeat_seq=' + repeat_seq + '\n']
                count += 1
                outfile.write('\t'.join(gff_line))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="trf.dat to GFF converter")
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