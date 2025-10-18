#!/usr/bin/env python3
from Bio import SeqIO
import os
import sys


#Takes a fasta file and a Hyde sample map and orders the sequences in the same order as the sample map, and writes output to phylip
infile = sys.argv[1] #the input fasta file
outfile = sys.argv[2] #the sorted output phylip file
sample_map = sys.argv[3] #the sample map (two columns, tab delimited, fasta sample name in the first column, population name (ignored) in the second)

n_seqs=0
seq_length=0

with open(infile, "r") as in_handle:
        for record in SeqIO.parse(in_handle, "fasta"):
                n_seqs+=1
                seq_length=len(record)
in_handle.close()


records = SeqIO.index(infile, "fasta")

fout = open(outfile,"w")
fout.write(" "+str(n_seqs)+"\t"+str(seq_length)+"\n")

with open(sample_map) as f:
    for line in f:
        clade=line.split()
        sample=clade[0]
        record=records[sample]
        fout.write(record.id+"\t"+str(record.seq)+"\n")
fout.close()
