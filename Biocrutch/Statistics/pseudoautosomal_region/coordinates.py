#!/usr/bin/env python3
__author__ = 'tomarovsky'
'''
Class for determining the coordinates of the pseudoautosomal region.
'''
from Biocrutch.Statistics.coverage_statistics.coverage_metrics import CoveragesMetrics
from collections import Counter


class Coordinates():
    def __init__(self, data, whole_genome_value, deviation_percent):
        self.data = data
        self.minimum_coverage = whole_genome_value - (whole_genome_value / 100 * deviation_percent)
        self.maximum_coverage = whole_genome_value + (whole_genome_value / 100 * deviation_percent)

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
        start_flag = True
        for lst in range(len(coordinates_list)):
            if not result:
                result.append(coordinates_list[lst])
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

    @staticmethod
    def median_coordinate_filter(coordinates: list, median_list: list, whole_genome_value: int, deviation_percent) -> list:
        minimum_coverage = whole_genome_value - (whole_genome_value / 100 * deviation_percent)
        # maximum_coverage = whole_genome_value + (whole_genome_value / 100 * deviation_percent)
        result = []
        empty_list = True
        median_index = -1
        
        for lst in range(len(coordinates)):
            if empty_list:
                result.append(coordinates[lst])
                empty_list = False
                continue
            median_index += 1
            if median_list[median_index] >= minimum_coverage:  # and coverage_value < maximum_coverage:
                result.append(coordinates[lst])
            continue
        return result


    def pseudocoordinates(self,
                          coverage_column_name: int,
                          window_column_name: int,
                          repeat_window_number: int) -> list:
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

            if coverage_value > self.minimum_coverage:  # and coverage_value < self.maximum_coverage:
                repeat_window += 1
                if repeat_window == repeat_window_number and start_coordinate is None:
                    start_coordinate = current_window - repeat_window + 1  # * args.window_size
                    repeat_window = 0
            elif start_coordinate is not None and (coverage_value <= self.minimum_coverage or coverage_value >= self.maximum_coverage):
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
        if between_regions_coverage_dict:
            median_between_regions.append(CoveragesMetrics(between_regions_coverage_dict).median_value())
            between_regions_coverage_dict.clear()
        
        print(coordinates)
        print(median_between_regions)
        print(self.minimum_coverage)
        
        # median concat
        #перепроверить
        result = []
        m = True
        for i in range(len(median_between_regions)):
            if median_between_regions[i] > self.minimum_coverage:  # and coverage_value < self.maximum_coverage:
                if m:
                    result.append(coordinates[i])
                if i != (len(median_between_regions)-1):
                    result.append(coordinates[i+1])
                    m = False
            else:
                m = True
                    
        return result