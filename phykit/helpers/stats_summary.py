import statistics as stat
import sys

import numpy as np


def calculate_summary_statistics_from_arr(
    arr    
):
    """
    calcuate summary statistics for an input list
    """
    stats = dict(
        mean               = stat.mean(arr),
        median             = stat.median(arr),
        twenty_fifth       = np.percentile(arr, 25),
        seventy_fifth      = np.percentile(arr, 75),
        minimum            = np.min(arr),
        maximum            = np.max(arr),
        standard_deviation = stat.stdev(arr),
        variance           = stat.variance(arr)
    )

    return stats

def calculate_summary_statistics_from_dict(
    dat: dict    
):
    """
    calcuate summary statistics for a dictionary
    """
    stats = dict(
        mean=stat.mean([*dat.values()]),
        median=stat.median([*dat.values()]),
        twenty_fifth=np.percentile([*dat.values()], 25),
        seventy_fifth=np.percentile([*dat.values()], 75),
        minimum=np.min([*dat.values()]),
        maximum=np.max([*dat.values()]),
        standard_deviation=stat.stdev([*dat.values()]),
        variance=stat.variance([*dat.values()])
    )

    return stats

def print_summary_statistics(
    stats: list
):
    """
    """
    try:
        print(f"mean: {round(stats['mean'], 4)}")
        print(f"median: {round(stats['median'], 4)}")
        print(f"25th percentile: {round(stats['twenty_fifth'], 4)}")
        print(f"75th percentile: {round(stats['seventy_fifth'], 4)}")
        print(f"minimum: {round(stats['minimum'], 4)}")
        print(f"maximum: {round(stats['maximum'], 4)}")
        print(f"standard deviation: {round(stats['standard_deviation'], 4)}")
        print(f"variance: {round(stats['variance'], 4)}")
    except BrokenPipeError:
        pass
