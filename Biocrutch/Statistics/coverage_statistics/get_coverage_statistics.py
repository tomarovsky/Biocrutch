#!/usr/bin/env python3
__author__ = 'tomarovsky'
from Biocrutch.Statistics.coverage_statistics.coverage_metrics import CoveragesMetrics
from collections import Counter
import pandas as pd


class GetCoverageStatistics:
    def __init__(self, data, output=False):
        self.data = data
        self.output = output

    def get_all_statistics(self, frame_size):
        # To calculate all statistics. The output is several files. (NO overlapping windows yet)
        # created dataframe for whole genome and scaffolds stats
        df_scaffolds = pd.DataFrame(columns=['median', 'average', 'max', 'min'])
        df_nonoverlapping_frames = pd.DataFrame(columns=['scaffold', 'frame', 'median', 'average', 'max', 'min'])
        df_whole_genome = pd.DataFrame(columns=['median', 'average', 'max', 'min'])

        genome_coverages_amounts_dict = Counter()
        genome_line_counter = 0
        frame_coverages_amounts_dict = Counter()
        frame_line_counter = 0
        frame_id = -1
        index = 0
        scaffold_coverages_dict = Counter()

        previous_scaffold_name = None

        for line in self.data:
            line = line.rstrip().split('\t')
            # for whole genome
            genome_line_counter += 1
            genome_coverages_amounts_dict[int(line[2])] += 1

            if previous_scaffold_name == line[0] or previous_scaffold_name == None:
                frame_line_counter += 1
                frame_coverages_amounts_dict[int(line[2])] += 1
            else:
                frame_id = -1
                frame_line_counter = 1
                frame_coverages_amounts_dict.clear()
                frame_coverages_amounts_dict[int(line[2])] += 1

            # for each scaffold
            if previous_scaffold_name != line[0] and previous_scaffold_name != None:
                metrics = CoveragesMetrics(scaffold_coverages_dict)
                print('scaffolds metrics is being processing')
                df_scaffolds.loc[previous_scaffold_name] = [metrics.median_value(),
                                                            metrics.average_value(),
                                                            metrics.max_coverage_value(),
                                                            metrics.min_coverage_value()]
                scaffold_coverages_dict.clear()

            # for window (non-overlapping)
            if frame_line_counter == frame_size:
                index += 1
                frame_id += 1
                metrics = CoveragesMetrics(frame_coverages_amounts_dict)
                print('non-overlapping windows metrics is being processing')
                df_nonoverlapping_frames.loc[index] = [previous_scaffold_name, frame_id, metrics.median_value(),
                                                                            metrics.average_value(),
                                                                            metrics.max_coverage_value(),
                                                                            metrics.min_coverage_value()]
                frame_coverages_amounts_dict.clear()
                frame_line_counter = 0

            scaffold_coverages_dict[int(line[2])] += 1
            previous_scaffold_name = line[0]

        # processing residual data after a cycle for each scaffold and whole genome
        metrics = CoveragesMetrics(scaffold_coverages_dict)
        print('scaffolds metrics is being processing')
        df_scaffolds.loc[previous_scaffold_name] = [metrics.median_value(),
                                                    metrics.average_value(),
                                                    metrics.max_coverage_value(),
                                                    metrics.min_coverage_value()]
        metrics = CoveragesMetrics(genome_coverages_amounts_dict)
        print('whole genome metrics is being processing')
        df_whole_genome.loc['whole_genome'] = [metrics.median_value(),
                                                    metrics.average_value(),
                                                    metrics.max_coverage_value(),
                                                    metrics.min_coverage_value()]

        #for print dataframe to terminal
        print(df_scaffolds)
        print(df_nonoverlapping_frames)
        print(df_whole_genome)

        # create a report.csv
        if self.output:
            df_scaffolds.rename_axis('scaffold').reset_index().to_csv(self.output + "_scaffolds_stats.csv", sep='\t', index = False)
            df_nonoverlapping_frames.to_csv(self.output + '_' + str(frame_size) + "_windows_stats.csv", sep='\t', index = False)
            df_whole_genome.rename_axis('genome').reset_index().to_csv(self.output + "_whole_genome_stats.csv", sep='\t', index = False)


    def get_whole_genome_stats(self):
        df_whole_genome = pd.DataFrame(columns=['median', 'average', 'max', 'min'])
        genome_coverages_amounts_dict = Counter()
        genome_line_counter = 0

        for line in self.data:
            line = line.rstrip().split('\t')
            genome_line_counter += 1
            genome_coverages_amounts_dict[int(line[2])] += 1

        # processing residual data after a cycle
        metrics = CoveragesMetrics(genome_coverages_amounts_dict)
        print('whole genome metrics is being processing')
        df_whole_genome.loc['whole_genome'] = [metrics.median_value(),
                                                    metrics.average_value(),
                                                    metrics.max_coverage_value(),
                                                    metrics.min_coverage_value()]
        # for print to terminal
        print(df_whole_genome)
        # create a report.csv
        if self.output:
            df_whole_genome.rename_axis('genome').reset_index().to_csv(self.output + "_whole_genome_stats.csv", encoding='utf-8', sep='\t', index = False)


    def get_scaffolds_stats(self):
        df_scaffolds = pd.DataFrame(columns=['median', 'average', 'max', 'min'])
        scaffold_coverages_dict = Counter()
        previous_scaffold_name = None

        for line in self.data:
            line = line.rstrip().split('\t')
            if previous_scaffold_name != line[0] and previous_scaffold_name != None:
                metrics = CoveragesMetrics(scaffold_coverages_dict)
                print('scaffolds metrics is being processing')
                df_scaffolds.loc[previous_scaffold_name] = [metrics.median_value(),
                                                            metrics.average_value(),
                                                            metrics.max_coverage_value(),
                                                            metrics.min_coverage_value()]
                scaffold_coverages_dict.clear()
            scaffold_coverages_dict[int(line[2])] += 1
            previous_scaffold_name = line[0]
        # processing residual data after a cycle
        metrics = CoveragesMetrics(scaffold_coverages_dict)
        print('scaffolds metrics is being processing')
        df_scaffolds.loc[previous_scaffold_name] = [metrics.median_value(),
                                                    metrics.average_value(),
                                                    metrics.max_coverage_value(),
                                                    metrics.min_coverage_value()]
        #for print dataframe to terminal
        print(df_scaffolds)
        # create a report.csv
        if self.output:
            df_scaffolds.rename_axis('scaffold').reset_index().to_csv(self.output + "_scaffolds_stats.csv", encoding='utf-8', sep='\t', index = False)


    def get_nonoverlapping_windows_stats(self, frame_size):
        df_nonoverlapping_frames = pd.DataFrame(columns=['scaffold', 'frame', 'median', 'average', 'max', 'min'])

        frame_coverages_amounts_dict = Counter()
        frame_line_counter = 0
        frame_id = -1
        index = 0
        previous_scaffold_name = None

        for line in self.data:
            line = line.rstrip().split('\t')
            if previous_scaffold_name == line[0] or previous_scaffold_name == None:
                frame_line_counter += 1
                frame_coverages_amounts_dict[int(line[2])] += 1
            else:
                frame_id = -1
                frame_line_counter = 1
                frame_coverages_amounts_dict.clear()
                frame_coverages_amounts_dict[int(line[2])] += 1
            # for window (non-overlapping)
            if frame_line_counter == frame_size:
                index += 1
                frame_id += 1
                metrics = CoveragesMetrics(frame_coverages_amounts_dict)
                print('non-overlapping windows metrics is being processing')
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
        if self.output:
            df_nonoverlapping_frames.to_csv(self.output + '_' + str(args.frame_size) + "_windows_stats.csv", encoding='utf-8', sep='\t', index = False)


    def get_overlapping_windows_stats(self): # in developing
        # df_frames = pd.DataFrame(columns=['scaffold', 'frame', 'median', 'average', 'max', 'min'])

        # overlapping_frame_id = 0.5
        # overlapping_frame_coverages_amounts_dict = Counter()
        # overlapping_frame_line_counter = 0
        # previous_scaffold_name = None

        # for line in args.input:
        #     genome_line_counter += 1
        #     line = line.rstrip().split('\t')

        #     if genome_line_counter > int(args.frame_size / 2):
        #         overlapping_frame_line_counter += 1
        #         overlapping_frame_coverages_amounts_dict[line[2]] += 1

        #     elif overlapping_frame_line_counter == args.frame_size:
        #         overlapping_frame_id += 1
        #         metrics = CoveragesMetrics(overlapping_frame_coverages_amounts_dict)
        #         df_frames.loc[overlapping_frame_id] = [previous_scaffold_name, metrics.median_value(),
        #                                                 metrics.average_value(),
        #                                                 metrics.max_coverage_value(),
        #                                                 metrics.min_coverage_value()]
        #         overlapping_frame_coverages_amounts_dict.clear()
        #         overlapping_frame_line_counter = 0

        # in developing :)
        pass


