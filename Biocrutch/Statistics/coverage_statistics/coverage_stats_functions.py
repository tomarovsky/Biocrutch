__author__ = 'tomarovsky'
from collections import Counter


def count_window_stats(coverages_amounts_dict: Counter) -> list:
    '''
    Function for calculating median, average, maximum and minimum coverage.
    Works with a Counter from collections packege.
    Output is a list with metrics.
    '''
    sum_of_coverages = sum(key*value for key, value in coverages_amounts_dict.items())
    sum_values_coverages = sum(coverages_amounts_dict.values())
    mean = round(sum_of_coverages/sum_values_coverages, 2)

    keys_coverages = sorted(coverages_amounts_dict.keys())
    half_sum_values_coverages = sum_values_coverages // 2
    count = 0

    if sum_values_coverages % 2 != 0:
        for key in keys_coverages:
            count += coverages_amounts_dict[key]
            if count >= half_sum_values_coverages:
                median = key
                metrics = [median, mean, keys_coverages[-1], keys_coverages[0]]
                print('window_stats: ', metrics)
                return metrics
    else:
        for i in range(len(keys_coverages)):
            count += coverages_amounts_dict[keys_coverages[i]]
            if count == half_sum_values_coverages:
                genome_median = (keys_coverages[i] + keys_coverages[i+1])/2
                metrics = [genome_median, mean, keys_coverages[-1], keys_coverages[0]]
                print('window_stats: ', metrics)
                return metrics
            elif count > half_sum_values_coverages:
                genome_median = keys_coverages[i]
                metrics = [genome_median, mean, keys_coverages[-1], keys_coverages[0]]
                print('window_stats: ', metrics)
                return metrics