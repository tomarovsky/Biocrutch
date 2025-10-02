#!/usr/bin/env python3
__author__ = 'tomarovsky'

import pandas as pd
import sys

if len(sys.argv) != 3:
    print("Usage: python  script.py  <input from bedtools intersect -wao>  <output_file>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

df = pd.read_csv(input_file, sep='\t', header=None)
df.columns = ['scaffold1', 'start', 'end', 'scaffold2', 'start2', 'end2', 'value']

result = df.groupby(['scaffold1', 'start', 'end'])['value'].sum().reset_index()

result.to_csv(output_file, sep='\t', index=False, header=False)
