__author__ = 'tomarovsky'
# Program descripton
DESCRIPTION = ('This program annotates sequence.fasta. '
               'Provides various features.')

# Text buffering
HELP_BUFFERING = ('Text buffering. Default = None')

# Minimum contig length
HELP_MIN_CONTIG = ('Cutoff allows you to set the minimum contig length.'
               'for statistical analysis. Default >= 0 and 500.'
               'The number of parameters is any.')

# N/L statistics
HELP_NL_STATISTICS = ('The parameter sets the values for N and L statistics.'
               'Default N50/L50 and N75/L75.'
               'The number of parameters is any.')

# GC content
HELP_GC_CONTENT = ('The parameter sets the values for GC content.'
               'Default parameters is 0 and 500.'
               'The number of parameters is any.')

# Input/Output
HELP_OUTPUT = ('Add output file prefix. The output file will be created in the current directory.')
HELP_INPUT = ('Input files. One or more files to run the program.'
              'The script also accepts files with the extension .gz, .bz2')

#Print to terminal
HELP_PRINT = ('If it is necessary to disable the output to the terminal, write "FALSE"'
              'The default is TRUE.')