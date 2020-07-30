from collections import OrderedDict
import pandas as pd
import bz2
import gzip

class Fasta_opener:
    def __init__(self, path):
        self.path = path
        self.lengths = {}

    @staticmethod
    def metaopen(filename, flags, buffering=None):
        if not isinstance(filename, str):
            return filename
        elif filename[-3:] == ".gz":
            return gzip.open(filename, flags)
        elif filename[-4:] == ".bz2":
            return bz2.open(filename, flags)
        else:
            if buffering is not None:
                return open(filename, flags, buffering=buffering)
            else:
                return open(filename, flags)

    def parse_sequences(self, buffering=None) -> dict:
        """
        Parses a text file of genome sequences into a dictionary.
        Arguments:
          buffering: buffer text value
        """
        print("parse_sequences started")
        data_dict = OrderedDict()
        self.lengths = {}
        header = None
        f = self.metaopen(self.path, 'rt', buffering)
        for line in f:
            if line.startswith('>'):
                header = line[1:].split(' ')[0]
                data_dict[header] = []
            else:
                data_dict[header].append(line[:-1])
        f.close()
        for name in data_dict:
            data_dict[name] = ''.join(data_dict[name])
            self.lengths[name] = len(data_dict[name])
        return data_dict

    def lengths_to_frame(self):
        """lengths to pandas dataframe"""
        print("lengths_to_frame started")
        lengths_df = pd.DataFrame.from_dict({'lengths' : self.lengths})
        return lengths_df


if __name__ == "__main__":
    data = Fasta_opener('data/sample.fasta.gz')
    print(data.parse_sequences(5000000))
    print(data.lengths_to_frame())
