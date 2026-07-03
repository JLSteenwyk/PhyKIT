from math import sqrt


class _LazyNumpy:
    _module = None

    def _load(self):
        if self._module is None:
            import numpy

            self._module = numpy
        return self._module

    def __getattr__(self, name):
        value = getattr(self._load(), name)
        setattr(self, name, value)
        return value


np = _LazyNumpy()

_NO_VALUES_MESSAGE = (
    "There are no values to calculate summary statistics for.\n\n"
    "Double check that the input alignment/phylogeny contains\n"
    "the properties you want to calculate summary statistics for."
)
_SMALL_SEQUENCE_STATS_MAX = 128
_SMALL_STATS_FALLBACK = object()


def _print_no_values_message():
    print(_NO_VALUES_MESSAGE)


def _no_values_result():
    _print_no_values_message()
    return None


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


def _linear_percentile(sorted_values, position):
    lower = int(position)
    upper = lower + 1
    if upper >= len(sorted_values):
        return sorted_values[lower]
    fraction = position - lower
    if fraction == 0:
        return sorted_values[lower]
    return sorted_values[lower] * (1.0 - fraction) + sorted_values[upper] * fraction


def _calculate_small_sequence_statistics(values):
    if not isinstance(values, (list, tuple)):
        return _SMALL_STATS_FALLBACK

    count = len(values)
    if count < 2 or count > _SMALL_SEQUENCE_STATS_MAX:
        return _SMALL_STATS_FALLBACK

    try:
        if count == 2:
            first_value, last_value = values
            all_integer = type(first_value) is int and type(last_value) is int
            if last_value < first_value:
                first_value, last_value = last_value, first_value

            if first_value == last_value:
                if all_integer:
                    mean = first_value
                    median = first_value
                    quartile = float(first_value)
                else:
                    mean = first_value
                    median = first_value
                    quartile = first_value
                return dict(
                    mean=mean,
                    median=median,
                    twenty_fifth=quartile,
                    seventy_fifth=quartile,
                    minimum=first_value,
                    maximum=last_value,
                    standard_deviation=0.0,
                    variance=0.0,
                )

            mean = (first_value + last_value) / 2
            median = mean
            span = last_value - first_value
            twenty_fifth = first_value + span * 0.25
            seventy_fifth = first_value + span * 0.75
            delta = first_value - mean
            variance = 2.0 * delta * delta
            if all_integer:
                mean = _integer_if_exact(mean)
                median = _integer_if_exact(median)
            return dict(
                mean=mean,
                median=median,
                twenty_fifth=twenty_fifth,
                seventy_fifth=seventy_fifth,
                minimum=first_value,
                maximum=last_value,
                standard_deviation=sqrt(variance),
                variance=variance,
            )

        sorted_values = sorted(values)
        first_value = sorted_values[0]
        last_value = sorted_values[-1]
        all_integer = type(first_value) is int and type(last_value) is int
        if all_integer:
            for index in range(1, count - 1):
                if type(sorted_values[index]) is not int:
                    all_integer = False
                    break
        if first_value == last_value:
            if all_integer:
                mean = first_value
                median = first_value
                quartile = float(first_value)
            else:
                mean = first_value
                median = first_value
                quartile = first_value
            return dict(
                mean=mean,
                median=median,
                twenty_fifth=quartile,
                seventy_fifth=quartile,
                minimum=first_value,
                maximum=last_value,
                standard_deviation=0.0,
                variance=0.0,
            )

        mean = sum(sorted_values) / count
        median = _linear_percentile(sorted_values, (count - 1) * 0.5)
        twenty_fifth = _linear_percentile(sorted_values, (count - 1) * 0.25)
        seventy_fifth = _linear_percentile(sorted_values, (count - 1) * 0.75)
        sum_squared_deviations = 0.0
        for value in sorted_values:
            delta = value - mean
            sum_squared_deviations += delta * delta
        variance = sum_squared_deviations / (count - 1)
    except (TypeError, ValueError):
        return _SMALL_STATS_FALLBACK

    if all_integer:
        mean = _integer_if_exact(mean)
        median = _integer_if_exact(median)

    return dict(
        mean=mean,
        median=median,
        twenty_fifth=twenty_fifth,
        seventy_fifth=seventy_fifth,
        minimum=first_value,
        maximum=last_value,
        standard_deviation=sqrt(variance),
        variance=variance,
    )


def _calculate_summary_statistics(values):
    if _has_too_few_values(values):
        return _no_values_result()

    small_stats = _calculate_small_sequence_statistics(values)
    if small_stats is not _SMALL_STATS_FALLBACK:
        return small_stats

    arr = np.asarray(values)
    if arr.size < 2:
        return _no_values_result()

    first_value = arr.flat[0]
    minimum = None
    maximum = None
    if first_value == arr.flat[-1]:
        minimum = arr.min()
        maximum = arr.max()
    if minimum is not None and minimum == maximum:
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
    if minimum is None:
        minimum = arr.min()
    if maximum is None:
        maximum = arr.max()
    variance = arr.var(ddof=1)
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
        return _no_values_result()

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
