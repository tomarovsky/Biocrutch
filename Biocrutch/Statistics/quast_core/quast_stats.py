__author__ = "tomarovsky"
from Biocrutch.Parsers.fasta_opener import Fasta_opener
import pandas as pd
import numpy as np
import re


class Quast_core:
    def __init__(self, sequences_dict, lengths_frame):
        self.sequences_dict = sequences_dict
        self.df = lengths_frame

    def contig_count(self, min_contig) -> int:
        print("contig_count started")
        return len(self.df[self.df["lengths"] >= min_contig].index)

    def largest_contig_lengh(self) -> int:
        print("largest_contig_length started")
        return self.df["lengths"].max()

    def total_length(self, min_contig) -> int:
        print("total_length started")
        return self.df[self.df["lengths"] >= min_contig]["lengths"].sum()

    def n_amount(self, min_contig):
        print("n_amount started")
        count = 0
        index_list = self.df[self.df["lengths"] >= min_contig].index
        for i in index_list:
            count += self.sequences_dict[i].upper().count("N")
        return count

    def gc_content(self, min_contig) -> float:
        print("gc_content started")
        count = 0
        index_list = self.df[self.df["lengths"] >= min_contig].index
        for i in index_list:
            count += self.sequences_dict[i].upper().count("G") + self.sequences_dict[i].upper().count("C")
        result = round(count / self.total_length(min_contig) * 100, 2)
        return result

    def n_l_statistics(self, percent: int, min_contig: int) -> list:
        print("n_l_statistics started")
        lengths = self.df[self.df["lengths"] >= min_contig]["lengths"].to_numpy(dtype=np.int64)
        if lengths.size == 0:
            return [None, None]

        lengths = np.sort(lengths)[::-1]
        total = lengths.sum()
        target = total * (percent / 100.0)

        csum = np.cumsum(lengths)
        idx = int(np.searchsorted(csum, target, side="left"))

        Nxx = int(lengths[idx])
        Lxx = idx + 1

        return [Nxx, Lxx]

    def n_l_contig_statistics(self, percent: int, min_contig: int) -> list:
        print("n_l_statistics contig started")
        contig_lengths = []

        for _, seq in self.sequences_dict.items():
            for match in re.finditer(r"[^Nn]+", seq):
                l = match.end() - match.start()
                if l >= min_contig:
                    contig_lengths.append(l)

        if not contig_lengths:
            return [None, None]

        lengths = np.sort(np.array(contig_lengths, dtype=np.int64))[::-1]

        total = lengths.sum()
        target = total * (percent / 100.0)

        csum = np.cumsum(lengths)
        idx = int(np.searchsorted(csum, target, side="left"))

        Nxx = int(lengths[idx])
        Lxx = idx + 1

        return [Nxx, Lxx]


if __name__ == "__main__":
    PATH = "data/sample.fna"
    fasta = Fasta_opener(PATH)
    sequences_dict = fasta.parse_sequences()
    metrics = Quast_core(sequences_dict, fasta.lengths_to_frame())
