__author__ = 'tomarovsky'
from Biocrutch.Routines.routine_functions import metaopen
from collections import OrderedDict
import pandas as pd


class Fasta_opener:
    def __init__(self, path):
        self.path = path
        self.lengths = {}

    def parse_sequences(self, buffering=None) -> dict:
        """
        Parses a text file of genome sequences into a dictionary.
        Arguments:
          buffering: buffer text value
        """
        print("parse_sequences started")
        data = OrderedDict()
        self.lengths = {}
        header = None
        f = metaopen(self.path, 'rt', buffering)
        for line in f:
            if line.startswith('>'):
                header = line[1:].split(' ')[0]
                data[header] = []
            else:
                data[header].append(line.rstrip())
        f.close()
        for name in data:
            data[name] = ''.join(data[name])
            self.lengths[name] = len(data[name])
        return data

    def parse_sequences_without_join(self, buffering=None) -> dict:
        print("parse_sequences started")
        data = OrderedDict()
        header = None
        f = metaopen(self.path, 'rt', buffering)
        for line in f:
            if line.startswith('>'):
                header = line
                data[header] = []
            else:
                data[header].append(line)
        f.close()
        for name in data:
            data[name] = ''.join(data[name])
        return data

    def lengths_to_frame(self):
        """lengths to pandas dataframe"""
        print("lengths_to_frame started")
        lengths_df = pd.DataFrame.from_dict({'lengths' : self.lengths})
        return lengths_df


if __name__ == "__main__":
    data = Fasta_opener('data/sample.fasta')
    print(data.parse_sequences(5000000))
    print(data.lengths_to_frame())
