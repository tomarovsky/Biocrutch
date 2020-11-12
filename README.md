# Biocrutch
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Сoncept scripts for bioinformatics.

### Content:

- [SRA_toolkit.py](https://github.com/etozhetoma/Biocrutch/blob/master/scripts/SRA_toolkit.py). The program parses the link from the Sequence Read Archive (SRA) and allows you to download reads in the sra format. The program also checks the integrity of the finished reads by parsing the required metrics of the source files.

- [quast_core.py](https://github.com/etozhetoma/Biocrutch/blob/master/scripts/quast_core.py). The program is an alternative to the publicly available Quast program. Its main differences are: 
    1. adding only the necessary cutoffs.
    2. counting missing N values.
    3. the output of the program is a pandas dataframe used for further analysis.
    4. output to a convenient csv file format.

- [get_coverage_statistics.py](https://github.com/etozhetoma/Biocrutch/blob/master/scripts/get_coverage_statistics.py). Script for calculating median, average, maximum and minimum coverage.
    1. calculate stats for whole genome
    2. calculate stats for each scaffold
    3. calculate stats in 100 kbp and 1 Mbp stacking windows

    Header-less tab-separated input file with 3 columns: scaffold_id, position(1-based), coverage.

- [get_pseudoautosomal_region.py](https://github.com/etozhetoma/Biocrutch/blob/master/scripts/get_pseudoautosomal_region.py). A script for determining the coordinates of the pseudo-autosomal region on the sex chromosome. The output is a BED file with the coordinates of the pseudoautosomal region.

- [ema_bin_to_fastq.py](https://github.com/etozhetoma/Biocrutch/blob/master/scripts/ema_bin_to_fastq.py). Combines Ema output files into reverse, forward and barcode-only file.

To use the Biocrutch package, you need to add the package path to PYTHONPATH.

Copyright (c) 2020 Andrey Tomarovsky
