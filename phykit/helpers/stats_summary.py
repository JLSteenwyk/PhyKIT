import statistics as stat


class _LazyNumpy:
    _module = None

    def _load(self):
        if self._module is None:
            import numpy

            self._module = numpy
        return self._module

    def __getattr__(self, name):
        return getattr(self._load(), name)


np = _LazyNumpy()

_NO_VALUES_MESSAGE = (
    "There are no values to calculate summary statistics for.\n\n"
    "Double check that the input alignment/phylogeny contains\n"
    "the properties you want to calculate summary statistics for."
)


def _print_no_values_message():
    print(_NO_VALUES_MESSAGE)


def _has_too_few_values(values):
    try:
        return len(values) < 2
    except TypeError:
        return False


def _python_scalar(value):
    return value.item() if hasattr(value, "item") else value


def _integer_if_exact(value):
    scalar = _python_scalar(value)
    if isinstance(scalar, float) and scalar.is_integer():
        return int(scalar)
    return scalar


def _calculate_summary_statistics(values):
    try:
        if _has_too_few_values(values):
            raise stat.StatisticsError

        arr = np.asarray(values)
        if arr.size < 2:
            raise stat.StatisticsError

        first_value = arr.flat[0]
        if first_value == arr.flat[-1] and np.all(arr == first_value):
            scalar = _python_scalar(first_value)
            if np.issubdtype(arr.dtype, np.integer):
                mean = _integer_if_exact(scalar)
                median = _integer_if_exact(scalar)
                quartile = float(scalar)
            else:
                mean = scalar
                median = scalar
                quartile = scalar
            return dict(
                mean=mean,
                median=median,
                twenty_fifth=quartile,
                seventy_fifth=quartile,
                minimum=scalar,
                maximum=scalar,
                standard_deviation=0.0,
                variance=0.0,
            )

        twenty_fifth, median, seventy_fifth = np.percentile(arr, [25, 50, 75])
        mean = arr.mean()
        minimum = np.min(arr)
        maximum = np.max(arr)
        variance = np.var(arr, ddof=1)
        standard_deviation = np.sqrt(variance)
        median = _python_scalar(median)
        if np.issubdtype(arr.dtype, np.integer):
            mean = _integer_if_exact(mean)
            median = _integer_if_exact(median)
        stats = dict(
            mean=_python_scalar(mean),
            median=median,
            twenty_fifth=_python_scalar(twenty_fifth),
            seventy_fifth=_python_scalar(seventy_fifth),
            minimum=_python_scalar(minimum),
            maximum=_python_scalar(maximum),
            standard_deviation=_python_scalar(standard_deviation),
            variance=_python_scalar(variance),
        )
    except stat.StatisticsError:
        _print_no_values_message()
        stats = None

    return stats


def calculate_summary_statistics_from_arr(arr):
    """
    calcuate summary statistics for an input list
    """
    return _calculate_summary_statistics(arr)


def calculate_summary_statistics_from_dict(dat: dict):
    """
    calcuate summary statistics for a dictionary
    """
    if len(dat) < 2:
        _print_no_values_message()
        return None

    try:
        arr = np.fromiter(dat.values(), dtype=float, count=len(dat))
    except (TypeError, ValueError):
        return _calculate_summary_statistics(list(dat.values()))
    return _calculate_summary_statistics(arr)


def print_summary_statistics(stats: list):
    """ """
    try:
        print(
            (
                "mean: %s\n"
                "median: %s\n"
                "25th percentile: %s\n"
                "75th percentile: %s\n"
                "minimum: %s\n"
                "maximum: %s\n"
                "standard deviation: %s\n"
                "variance: %s"
            )
            % (
                round(stats["mean"], 4),
                round(stats["median"], 4),
                round(stats["twenty_fifth"], 4),
                round(stats["seventy_fifth"], 4),
                round(stats["minimum"], 4),
                round(stats["maximum"], 4),
                round(stats["standard_deviation"], 4),
                round(stats["variance"], 4),
            )
        )
    except BrokenPipeError:
        pass
