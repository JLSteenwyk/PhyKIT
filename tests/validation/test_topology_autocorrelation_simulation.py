"""Independent expected-statistic checks and small reproducible spatial experiments."""

import importlib.util
import itertools
from pathlib import Path

import numpy as np
import pytest

from phykit.helpers.topology_autocorrelation import analyze, read_classified

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("autocorrelation_calibration",
                                            ROOT / "benchmarks/calibrate_topology_autocorrelation.py")
CALIBRATION = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CALIBRATION)


@pytest.mark.validation
@pytest.mark.parametrize("length", [0, 250])
def test_analytic_target_matches_independent_small_pair_enumeration(length):
    x = np.array([0, 100, 300, 500, 1000])
    p = np.tile([0.5, 0.3, 0.2], (5, 1))
    if not length:
        p[0], p[4] = [0.8, 0.1, 0.1], [0.1, 0.1, 0.8]
    mask = np.array([True, True, False, True, True])
    expected_joint, pair_distances = [], []
    for i, j in itertools.combinations(np.flatnonzero(mask), 2):
        distance = x[j] - x[i]
        joint = p[i] * p[j]
        if length:
            joint += p[i] * (1 - p[i]) * np.exp(-distance / length)
        expected_joint.append(joint)
        pair_distances.append(distance)
    expected_joint = np.array(expected_joint)
    expected = []
    for lower, upper in itertools.pairwise(CALIBRATION.EDGES):
        selected = (np.array(pair_distances) >= lower) & (np.array(pair_distances) < upper)
        specific = expected_joint[selected].mean(axis=0) - expected_joint.mean(axis=0)
        expected.append([specific.sum(), *specific])
    np.testing.assert_allclose(CALIBRATION.target([(x, p, mask)], length), expected, atol=1e-14)


@pytest.mark.validation
def test_equal_totals_different_spatial_patterns():
    fixture = ROOT / "tests/sample_files/topology_autocorrelation"
    a, _ = analyze(read_classified(fixture / "clustered.tsv"), [0, 101])
    b, _ = analyze(read_classified(fixture / "interleaved.tsv"), [0, 101])
    assert a["reference_estimates"][0]["excess"] == pytest.approx(6 / 11)
    assert b["reference_estimates"][0]["excess"] == pytest.approx(-3 / 11)
    assert a["chromosome_diagnostics"] == b["chromosome_diagnostics"]


@pytest.mark.validation
def test_simulation_estimator_matches_iid_and_clustered_targets():
    for name, length in [("iid_irregular_missing", 0), ("clustered_irregular_missing", 250)]:
        chromosomes = CALIBRATION.geometry(name, 1833)
        target = CALIBRATION.target(chromosomes, length)
        results = []
        for seed in range(20):
            rows = CALIBRATION.simulate(chromosomes, length, np.random.default_rng(seed))
            estimate, _ = analyze(rows, CALIBRATION.EDGES)
            results.append(np.array([r["excess"] for r in estimate["reference_estimates"]]).reshape(2, 4))
        results = np.array(results)
        error = np.abs(results.mean(axis=0) - target)
        # Four Monte Carlo standard errors plus a small absolute rounding margin.
        assert np.all(error < 4 * results.std(axis=0, ddof=1) / np.sqrt(len(results)) + 0.001)
        if length:
            assert results.mean(axis=0)[0, 0] > results.mean(axis=0)[1, 0] > 0


@pytest.mark.validation
def test_independent_nonstationary_labels_produce_confounded_excess():
    for name in ("frequency_gradient", "spatial_missingness"):
        target = CALIBRATION.target(CALIBRATION.geometry(name, 991), 0)
        assert target[0, 0] > 0.04
