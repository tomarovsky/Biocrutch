#!/usr/bin/env python3
__author__ = 'tomarovsky'
from Biocrutch.Statistics.coverage_statistics.CoverageMetrics import CoveragesMetrics
from collections import Counter
import pandas as pd


class GenomecovCoverageStatistics:
    def __init__(self, data, output, tool_name):
        self.data = data
        self.output = output # if not output.endswith(".csv") else output[-4:]
        self.tool_name = tool_name
        

    def get_whole_genome_stats(self):
        df_whole_genome = pd.DataFrame(columns=['median', 'average', 'max', 'min'])
        genome_coverages_amounts_dict = Counter()
        genome_line_counter = 0

        for line in self.data:
            line = line.rstrip().split('\t')
            genome_line_counter += 1
            genome_coverages_amounts_dict[float(line[2])] += 1

        # processing residual data after a cycle
        metrics = CoveragesMetrics(genome_coverages_amounts_dict)
        # print('whole genome metrics is being processing')
        df_whole_genome.loc['whole_genome'] = [metrics.median_value(),
                                               metrics.average_value(),
                                               metrics.max_coverage_value(),
                                               metrics.min_coverage_value()]
        # for print to terminal
        print(df_whole_genome)
        # create a report.csv
        df_whole_genome.rename_axis('genome').reset_index().to_csv(self.output + "_whole_genome_stats.csv",
                                                                   encoding='utf-8', sep='\t', index = False)
    
    
    def get_scaffolds_stats(self):
        df_scaffolds = pd.DataFrame(columns=['median', 'average', 'max', 'min'])
        scaffold_coverages_dict = Counter()
        previous_scaffold_name = None

        for line in self.data:
            line = line.rstrip().split('\t')
            if previous_scaffold_name != line[0] and previous_scaffold_name != None:
                metrics = CoveragesMetrics(scaffold_coverages_dict)
                # print('scaffolds metrics is being processing')
                df_scaffolds.loc[previous_scaffold_name] = [metrics.median_value(),
                                                            metrics.average_value(),
                                                            metrics.max_coverage_value(),
                                                            metrics.min_coverage_value()]
                scaffold_coverages_dict.clear()
            scaffold_coverages_dict[float(line[2])] += 1
            previous_scaffold_name = line[0]
        # processing residual data after a cycle
        metrics = CoveragesMetrics(scaffold_coverages_dict)
        # print('scaffolds metrics is being processing')
        df_scaffolds.loc[previous_scaffold_name] = [metrics.median_value(),
                                                    metrics.average_value(),
                                                    metrics.max_coverage_value(),
                                                    metrics.min_coverage_value()]
        #for print dataframe to terminal
        print(df_scaffolds)
        # create a report.csv
        df_scaffolds.rename_axis('scaffold').reset_index().to_csv(self.output + "_scaffolds_stats.csv",
                                                                  encoding='utf-8', sep='\t', index = False)


    def get_nonoverlapping_windows_stats(self, frame_size):
        df_nonoverlapping_frames = pd.DataFrame(columns=['scaffold', 'frame', 'median', 'average', 'max', 'min'])

        frame_coverages_amounts_dict = Counter()
        frame_line_counter = 0
        frame_id = -1
        index = 0
        previous_scaffold_name = None

        for line in self.data:
            line = line.rstrip().split('\t')
            if previous_scaffold_name == line[0] or previous_scaffold_name is None:
                frame_line_counter += 1
                frame_coverages_amounts_dict[float(line[2])] += 1
            else:
                frame_id = -1
                frame_line_counter = 1
                frame_coverages_amounts_dict.clear()
                frame_coverages_amounts_dict[float(line[2])] += 1
            # for window (non-overlapping)
            if frame_line_counter == frame_size:
                index += 1
                frame_id += 1
                metrics = CoveragesMetrics(frame_coverages_amounts_dict)
                # print('non-overlapping windows metrics is being processing')
                df_nonoverlapping_frames.loc[index] = [previous_scaffold_name, frame_id, metrics.median_value(),
                                                                            metrics.average_value(),
                                                                            metrics.max_coverage_value(),
                                                                            metrics.min_coverage_value()]
                frame_coverages_amounts_dict.clear()
                frame_line_counter = 0
            previous_scaffold_name = line[0]

        #for print dataframe to terminal
        print(df_nonoverlapping_frames)
        # create a report.csv
        df_nonoverlapping_frames.to_csv(self.output + '_' + str(frame_size) + "_windows_stats.csv",
                                        encoding='utf-8', sep='\t', index = False)

    def get_universal_windows_stats(self, frame_size, frame_shift):
        df_overlapping_frames = pd.DataFrame(columns=['scaffold', 'frame', 'median', 'average', 'max', 'min'])
        data = self.data.readlines()
        coverages_dict = Counter()
        frame_id = -1 # for numbering from 0
        index = 0
        gap_counter = 0

        for ln in range(0, len(data), frame_shift):
            try:
                scaffold_name = data[ln + gap_counter].rstrip().split('\t')[0]
                last_scaffold_name = data[ln - 1 + frame_size + gap_counter].rstrip().split('\t')[0]
                # next_scaffold_name = data[ln + frame_size + gap_counter].rstrip().split('\t')[0]
                # print("actual scaffold:", scaffold_name)
                # print("last scaffold:", last_scaffold_name)

                if scaffold_name != last_scaffold_name:
                    for gap in range(frame_size + 1):
                        scaffold_name = data[ln + gap_counter + gap].rstrip().split('\t')[0]
                        # print("ELSE: now and next", scaffold_name, last_scaffold_name)
                        if scaffold_name == last_scaffold_name:
                            gap_counter += gap
                            # print("gap_counter:", gap_counter)
                            scaffold_name = data[ln + gap_counter].rstrip().split('\t')[0]
                            last_scaffold_name = data[ln - 1 + frame_size + gap_counter].rstrip().split('\t')[0]
                            # next_scaffold_name = data[ln + frame_size + gap_counter].rstrip().split('\t')[0]
                            break

                if scaffold_name == last_scaffold_name:
                    for j in range(frame_size):
                        line = data[ln + gap_counter + j].rstrip().split('\t')
                        # print("data line:", line)
                        coverages_dict[float(line[2])] += 1
                        if j == frame_size - 1:
                            index += 1
                            frame_id += 1
                            metrics = CoveragesMetrics(coverages_dict)
                            # print('universal windows metrics is being processing')
                            df_overlapping_frames.loc[index] = [scaffold_name,
                                                                frame_id, 
                                                                metrics.median_value(),
                                                                metrics.average_value(),
                                                                metrics.max_coverage_value(),
                                                                metrics.min_coverage_value()]
                            coverages_dict.clear()
            except IndexError:
                break
            
        #for print dataframe to terminal
        print(df_overlapping_frames)
        # create a report.csv
        df_overlapping_frames.to_csv(self.output + '_' + str(frame_size) + "_windows_stats.csv", encoding='utf-8', sep='\t', index = False)
