#!/usr/bin/env python3
__author__ = 'tomarovsky'
from collections import Counter


class CoveragesMetrics:
    '''
    Class for calculating median, average, maximum and minimum coverage.
    Works with a Counter from collections packege.
    '''
    def __init__(self, coverages_amounts_dict: Counter):
        self.coverages_amounts_dict = coverages_amounts_dict
        self.sum_values_coverages = sum(coverages_amounts_dict.values())

    def max_coverage_value(self) -> int:
        return max(self.coverages_amounts_dict.keys())

    def min_coverage_value(self) -> int:
        return min(self.coverages_amounts_dict.keys())

    def average_value(self) -> float:
        sum_of_coverages = sum(key*value for key, value in self.coverages_amounts_dict.items())
        average = sum_of_coverages/self.sum_values_coverages
        return average

    def median_value(self) -> int:
        keys_coverages = sorted(self.coverages_amounts_dict.keys())
        half_sum_values_coverages = self.sum_values_coverages // 2
        count = 0
        if self.sum_values_coverages % 2 != 0:
            for key in keys_coverages:
                count += self.coverages_amounts_dict[key]
                if count > half_sum_values_coverages:
                    median = key
                    return median
        else:
            for i in range(len(keys_coverages)):
                count += self.coverages_amounts_dict[keys_coverages[i]]
                if count == half_sum_values_coverages:
                    median = (keys_coverages[i] + keys_coverages[i+1]) / 2
                    return median
                elif count > half_sum_values_coverages:
                    median = keys_coverages[i]
                    return median

# bad draft function
# def count_window_stats(coverages_amounts_dict: Counter) -> list:
#     '''
#     Function for calculating median, average, maximum and minimum coverage.
#     Works with a Counter from collections packege.
#     Output is a list with metrics.
#     '''
#     sum_of_coverages = sum(key*value for key, value in coverages_amounts_dict.items())
#     sum_values_coverages = sum(coverages_amounts_dict.values())
#     mean = round(sum_of_coverages/sum_values_coverages, 2)

#     keys_coverages = sorted(coverages_amounts_dict.keys())
#     half_sum_values_coverages = sum_values_coverages // 2
#     count = 0

#     if sum_values_coverages % 2 != 0:
#         for key in keys_coverages:
#             count += coverages_amounts_dict[key]
#             if count > half_sum_values_coverages:
#                 median = key
#                 metrics = [median, mean, keys_coverages[-1], keys_coverages[0]]
#                 print('window_stats: ', metrics)
#                 return metrics
#     else:
#         for i in range(len(keys_coverages)):
#             count += coverages_amounts_dict[keys_coverages[i]]
#             if count == half_sum_values_coverages:
#                 median = (keys_coverages[i] + keys_coverages[i+1]) / 2
#                 metrics = [median, mean, keys_coverages[-1], keys_coverages[0]]
#                 print('window_stats: ', metrics)
#                 return metrics
#             elif count > half_sum_values_coverages:
#                 median = keys_coverages[i]
#                 metrics = [median, mean, keys_coverages[-1], keys_coverages[0]]
#                 print('window_stats: ', metrics)
#                 return metrics