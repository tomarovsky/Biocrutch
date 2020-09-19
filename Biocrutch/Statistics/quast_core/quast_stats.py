from Biocrutch.Parsers.fasta_opener import Fasta_opener
import pandas as pd

class Quast_core:
    def __init__(self, sequences_dict, lengths_frame):
        self.sequences_dict = sequences_dict
        self.df = lengths_frame

    def contig_count(self, min_contig=500) -> int:
        print("contig_count started")
        return len(self.df[self.df['lengths'] > min_contig].index)

    def largest_contig_lengh(self) -> int:
        print("largest_contig_lengh started")
        return self.df['lengths'].max()

    def total_length(self, min_contig=500) -> int:
        print("total_length started")
        return self.df[self.df['lengths'] > min_contig]['lengths'].sum()

    def n_amount(self, min_contig=500):
        print("n_amount started")
        count = 0
        index_list = self.df[self.df['lengths'] >= min_contig].index
        for i in index_list:
            count += self.sequences_dict[i].upper().count('N')
        return count

    def gc_content(self, min_contig=500) -> float:
        print("gc_content started")
        count = 0
        index_list = self.df[self.df['lengths'] >= min_contig].index
        for i in index_list:
            count += self.sequences_dict[i].upper().count('G') + self.sequences_dict[i].upper().count('C')
        result = round(count / self.total_length(min_contig) * 100, 2)
        return result

    def n_l_statistics(self, percent=50, min_contig=500) -> list:
        print("n_l_statistics started")
        mean = self.total_length(min_contig) // 100 * percent
        lengths = self.df[self.df['lengths'] > min_contig]['lengths'].sort_values(ascending=False)
        if lengths.empty:
            return [None, None]
        l_count = 0
        count = 0
        for l in lengths:
            count += l
            if count <= mean:
                l_count += 1
            else:
                l_count += 1
                n_count = l
                return [n_count, l_count]

if __name__ == "__main__":
    PATH = 'data/sample.fna'
    fasta = Fasta_opener(PATH)
    sequences_dict = fasta.parse_sequences()
    metrics = Quast_core(sequences_dict, fasta.lengths_to_frame())
    print(metrics.n_amount())
