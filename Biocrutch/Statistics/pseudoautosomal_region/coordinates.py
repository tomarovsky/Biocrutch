#!/usr/bin/env python3
__author__ = 'tomarovsky'
'''
Class for determining the coordinates of the pseudoautosomal region.
'''
from Biocrutch.Statistics.coverage_statistics.coverage_metrics import CoveragesMetrics
from collections import Counter


class Coordinates():
    def __init__(self, data):
        self.data = data

    def pseudocoordinates(self,
                          whole_genome_value: int,
                          deviation_percent: int,
                          coverage_column_name: int,
                          window_column_name: int,
                          repeat_window_number: int) -> list:
        deviation = whole_genome_value / 100 * deviation_percent
        minimum_coverage = whole_genome_value - deviation  # 25.5
        maximum_coverage = whole_genome_value + deviation  # 42.5

        coordinates = []
        repeat_window = 0
        start_coordinate = None

        between_regions_coverage_dict = Counter()
        median_between_regions = []
        between_region_flag = None

        for ln in self.data:
            line = ln.rstrip().split("\t")
            coverage_value = float(line[coverage_column_name])
            current_window = int(line[window_column_name])

            if between_region_flag:
                between_regions_coverage_dict[coverage_value] += 1

            if coverage_value > minimum_coverage:  # and coverage_value < maximum_coverage:
                repeat_window += 1
                if repeat_window == repeat_window_number and start_coordinate is None:
                    start_coordinate = current_window - repeat_window + 1  # * args.window_size
                    repeat_window = 0
            elif start_coordinate is not None and (coverage_value <= minimum_coverage or coverage_value >= maximum_coverage):
                stop_coordinate = current_window - 1  # * args.window_size
                coordinates.append([start_coordinate, stop_coordinate])
                if between_region_flag:
                    median_between_regions.append(CoveragesMetrics(between_regions_coverage_dict).median_value())
                    between_regions_coverage_dict.clear()
                if coordinates:
                    between_region_flag = True
                start_coordinate = None
                repeat_window = 0
            else:
                repeat_window = 0
        print(coordinates)
        print(median_between_regions)

        return coordinates

    @staticmethod
    def coordinates_between_regions(coordinates: list, between_regions_list: list) -> list:
        # input list [[start, stop], [start, stop]]
        last_element_index = len(coordinates) - 1
        for lst in range(len(coordinates)):
            start = coordinates[lst][-1]
            if lst == last_element_index:
                break
            stop = coordinates[lst + 1][0]
            between_regions_list.append([start, stop])
        return between_regions_list

    @staticmethod
    def distanse_coordinate_filter(coordinates_list: list, min_region_length: int) -> list:
        # input list [[start, stop], [start, stop]]
        result = []
        empty_list = True
        start_flag = True
        for lst in range(len(coordinates_list)):
            if empty_list:
                result.append(coordinates_list[lst])
                empty_list = False
                continue
            d_first = coordinates_list[lst - 1][1]
            d_second = coordinates_list[lst][0]
            distanse = d_second - d_first
            if distanse > min_region_length:
                result.append(coordinates_list[lst])
            else:
                if start_flag == True:
                    start = coordinates_list[lst - 1][0]
                    start_flag = False
                stop = coordinates_list[lst][1]
                if result[-1][1] < start:
                    result.append([start, stop])
                    start_flag = True
                else:
                    del result[-1][-1]
                    result[-1].append(stop)
        return result
