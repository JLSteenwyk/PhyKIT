import statistics as stat

import numpy as np


def calculate_summary_statistics(
    arr    
):
    """
    calcuate summary statistics for an input list
    """
    # TODO: break out stats calculations to a base function
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

def print_summary_statistics(
    stats: list
):
    """
    """
    print(f"mean: {stats['mean']}")
    print(f"median: {stats['median']}")
    print(f"25th percentile: {stats['twenty_fifth']}")
    print(f"75th percentile: {stats['seventy_fifth']}")
    print(f"minimum: {stats['minimum']}")
    print(f"maximum: {stats['maximum']}")
    print(f"standard deviation: {stats['standard_deviation']}")
    print(f"variance: {stats['variance']}")