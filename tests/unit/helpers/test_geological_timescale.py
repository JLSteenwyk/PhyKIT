import subprocess
import sys

import phykit.helpers.geological_timescale as timescale
from phykit.helpers.geological_timescale import get_timescale_for_range


def test_module_import_has_no_heavy_dependencies():
    code = """
import sys
import phykit.helpers.geological_timescale
assert "Bio" not in sys.modules
assert "numpy" not in sys.modules
assert "matplotlib" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_get_timescale_for_range_auto_epoch_matches_explicit_epoch():
    auto_intervals, auto_colors = get_timescale_for_range(40, "auto")
    epoch_intervals, epoch_colors = get_timescale_for_range(40, "epoch")

    assert auto_intervals == epoch_intervals
    assert auto_colors is epoch_colors is timescale.EPOCH_COLORS
    assert auto_intervals[0][0] == "Holocene"


def test_get_timescale_for_range_epoch_scans_table_once(monkeypatch):
    class CountingEpochs(list):
        iterations = 0

        def __iter__(self):
            self.iterations += 1
            return super().__iter__()

    epochs = CountingEpochs(timescale.EPOCHS)
    monkeypatch.setattr(timescale, "EPOCHS", epochs)

    intervals, colors = get_timescale_for_range(40, "epoch")

    assert intervals
    assert colors is timescale.EPOCH_COLORS
    assert epochs.iterations == 1


def test_get_timescale_for_range_epoch_stops_after_first_out_of_range_interval(
    monkeypatch,
):
    class EarlyStopEpochs:
        def __iter__(self):
            yield ("In range", 0.5, 0)
            yield ("Out of range", 2.0, 1.5)
            raise AssertionError("later intervals should not be scanned")

    monkeypatch.setattr(timescale, "EPOCHS", EarlyStopEpochs())

    intervals, colors = get_timescale_for_range(1, "epoch")

    assert intervals == [("In range", 0.5, 0)]
    assert colors is timescale.EPOCH_COLORS


def test_get_timescale_for_range_period_stops_after_first_out_of_range_interval(
    monkeypatch,
):
    class EarlyStopPeriods:
        def __iter__(self):
            yield ("In range", 0.5, 0)
            yield ("Out of range", 2.0, 1.5)
            raise AssertionError("later intervals should not be scanned")

    monkeypatch.setattr(timescale, "PERIODS", EarlyStopPeriods())

    intervals, colors = get_timescale_for_range(1, "period")

    assert intervals == [("In range", 0.5, 0)]
    assert colors is timescale.PERIOD_COLORS
