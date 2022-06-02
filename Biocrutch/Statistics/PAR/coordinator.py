#!/usr/bin/env python3
__author__ = 'tomarovsky'
'''
Class for determining the coordinates of the pseudoautosomal region.
'''
from Biocrutch.Statistics.coverage_statistics.CoverageMetrics import CoveragesMetrics
from collections import Counter


class Coordinator:
    def __init__(self, data, whole_genome_value, region_gap_size, deviation_percent):
        self.data = data
        self.whole_genome_value = whole_genome_value
        self.region_gap_size = region_gap_size
        self.minimum_coverage = whole_genome_value - (whole_genome_value / 100 * deviation_percent)
        self.maximum_coverage = whole_genome_value + (whole_genome_value / 100 * deviation_percent)

    def get_coordinates(self, window_size,
                    coverage_column_name: int,
                    window_column_name: int,
                    repeat_window_number: int) -> list:

        coordinates = []
        coverages_between_regions = []
        between_regions_flag = False
        median_between_regions_list = []
        repeat_window = 0
        start_coordinate = None

        for line in self.data:
            line = line.rstrip().split("\t")
            coverage_value = float(line[coverage_column_name])
            current_window = int(line[window_column_name])

            if between_regions_flag:
                coverages_between_regions.append(coverage_value)

            if coverage_value >= self.minimum_coverage:
                repeat_window += 1
                if repeat_window == repeat_window_number and start_coordinate is None:
                    start_coordinate = (current_window - repeat_window + 1) * window_size
                    if between_regions_flag:
                        coverages_between_regions = coverages_between_regions[:-repeat_window]
                        if len(coverages_between_regions) >= self.region_gap_size:
                            # the median of the section between regions, which is less than region_gap_size, is considered acceptable
                            between_regions_coverage_dict = Counter()
                            for i in coverages_between_regions:
                                between_regions_coverage_dict[i] += 1
                            median_between_regions_list.append(CoveragesMetrics(between_regions_coverage_dict).median_value())
                            between_regions_coverage_dict.clear()
                            coverages_between_regions = []
                        else:
                            # if the distance between regions < minimum_coverage, then we consider the admissible median.
                            median_between_regions_list.append(self.minimum_coverage)
                        between_regions_flag = False
            elif coverage_value < self.minimum_coverage and start_coordinate is not None:
                stop_coordinate = current_window * window_size
                coordinates.append([start_coordinate, stop_coordinate])
                between_regions_flag = True
                start_coordinate = None
                repeat_window = 0
            else:
                repeat_window = 0

        if start_coordinate is not None:
            stop_coordinate = current_window * window_size
            coordinates.append([start_coordinate, stop_coordinate])

        return coordinates, median_between_regions_list
