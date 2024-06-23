# Biocrutch
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Various bioinformatics scripts.

### Content:

- [SRAtoolkit](https://github.com/tomarovsky/Biocrutch/blob/master/scripts/Auto/SRA_toolkit.py). The program parses the link from the Sequence Read Archive (SRA) and allows you to download reads in the sra format. The program also checks the integrity of the finished reads by parsing the required metrics of the source files.
- [Pseudoautosomal region](https://github.com/tomarovsky/Biocrutch/blob/master/scripts/PAR/pseudoautosomal_region.py). A script for determining the coordinates of the pseudo-autosomal region on the sex chromosome. The output is a BED file with the coordinates of the pseudoautosomal region.
- [RepeatMasking scripts](https://github.com/tomarovsky/Biocrutch/tree/master/scripts/RepeatMasking). Scripts for converting TRF, RepeatMasker and WindowMasker output to GFF format.
- [EMA to FASTQ](https://github.com/tomarovsky/Biocrutch/blob/master/scripts/10x/ema_bin_to_fastq.py). Combines Ema output BINfiles into reverse FASTQ, forward FASTQ and barcode-only files.
- [QuastCore](https://github.com/tomarovsky/Biocrutch/blob/master/scripts/Statistics/quast_core.py). The program is an alternative to the publicly available Quast program. Its main differences are:
    1. adding only the necessary cutoffs.
    2. counting missing N values.
    3. the output of the program is a pandas dataframe used for further analysis.
    4. output to a convenient csv file format. And others.
- [Coverage statistics](https://github.com/tomarovsky/Biocrutch/blob/master/scripts/Coverage/coverage_statistics.py). Script for calculating median, average, maximum and minimum coverage. Script works with the output of Bedtools Genomecov and Mosdepth programs.
    1. calculate stats for whole genome.
    2. calculate stats for each scaffold.
    3. calculate stats stacking windows.
- [PSMC data combine](https://github.com/tomarovsky/Biocrutch/blob/master/scripts/Auto/psmc_data_combine.py). The script combines data from several [PSMC](https://github.com/lh3/psmc) outputs to draw multiple demographic population histories on a single graph.

and others...

To use the Biocrutch package, you need to add the package path to PYTHONPATH.

Copyright (c) 2020 Andrey Tomarovsky
