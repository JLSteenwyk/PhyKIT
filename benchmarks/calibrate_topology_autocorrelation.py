"""Reproducible coverage experiments with independently derived spatial targets."""

import argparse
import itertools
import json
import platform
import time
from pathlib import Path

import numpy as np
from scipy.stats import binomtest

from phykit.helpers.topology_autocorrelation import METRICS, analyze, uncertainty
from phykit.helpers.topology_landscape import CLASSES

EDGES = [0, 201, 501]
BLOCKS = [5000, 10000]
SCENARIOS = ("iid_regular", "iid_irregular_missing", "iid_imbalanced",
             "clustered_regular", "clustered_irregular_missing", "frequency_gradient",
             "spatial_missingness")


def geometry(name, seed):
    rng = np.random.default_rng(seed)
    chromosomes = []
    for c, n in enumerate((6000, 4000)):
        irregular = "irregular" in name or name in ("frequency_gradient", "spatial_missingness")
        x = np.cumsum(rng.integers(25, 176, n) if irregular else np.full(n, 100))
        p = np.asarray(([0.5, 0.3, 0.2], [0.2, 0.3, 0.5])[c])
        if name == "iid_regular":
            p = np.full(3, 1 / 3)
        if name == "iid_imbalanced":
            p = np.asarray([0.8, 0.15, 0.05])
        probabilities = np.tile(p, (n, 1))
        if name in ("frequency_gradient", "spatial_missingness"):
            t = np.linspace(0, 1, n)
            probabilities = np.column_stack((0.8 - 0.6 * t, np.full(n, 0.15), 0.05 + 0.6 * t))
        observed = rng.random(n) >= (0.2 if "missing" in name else 0)
        if name == "spatial_missingness":
            observed = rng.random(n) >= np.where(x < np.median(x), 0.1, 0.7)
        chromosomes.append((x, probabilities, observed))
    return chromosomes


def target(chromosomes, length):
    """Expected finite-statistic excess conditional on positions and missing mask.

    Under a spatial reset process, P(Y_i=Y_j=k) = p_k^2 +
    p_k(1-p_k) exp(-distance/length). The global baseline expectation uses
    every unordered pair, not the marginal-frequency limit sum p_k^2.
    """
    numerator = np.zeros((len(EDGES) - 1, 3))
    denominator = np.zeros(len(EDGES) - 1)
    for x, probabilities, observed in chromosomes:
        x, probabilities = x[observed], probabilities[observed]
        n = len(x)
        global_joint = (probabilities.sum(axis=0)**2 - (probabilities**2).sum(axis=0)) / (n * (n - 1))
        if length:
            previous = total = 0.
            for distance in np.diff(x):
                previous = (previous + 1) * np.exp(-distance / length)
                total += previous
            p = probabilities[0]
            global_joint += p * (1 - p) * total / (n * (n - 1) / 2)
        for i, position in enumerate(x):
            for b, (lower, upper) in enumerate(itertools.pairwise(EDGES)):
                first = max(i + 1, int(np.searchsorted(x, position + lower)))
                last = int(np.searchsorted(x, position + upper))
                if last <= first:
                    continue
                joint = probabilities[i] * probabilities[first:last]
                if length:
                    p = probabilities[i]
                    joint += p * (1 - p) * np.exp(-(x[first:last] - position) / length)[:, None]
                numerator[b] += joint.sum(axis=0) - (last - first) * global_joint
                denominator[b] += last - first
    specific = numerator / denominator[:, None]
    return np.column_stack((specific.sum(axis=1), specific))


def simulate(chromosomes, length, rng):
    rows = []
    for c, (x, p, observed) in enumerate(chromosomes):
        labels = (rng.random(len(x))[:, None] > np.cumsum(p, axis=1)).sum(axis=1)
        if length:
            retained = rng.random(len(x) - 1) < np.exp(-np.diff(x) / length)
            for i in range(1, len(x)):
                if retained[i - 1]:
                    labels[i] = labels[i - 1]
        for i, (position, label, keep) in enumerate(zip(x, labels, observed)):
            rows.append({"reference": "ref", "chromosome": f"chr{c + 1}", "gene_id": f"c{c}g{i}",
                         "start": int(position), "end": int(position) + 1, "anchor": int(position),
                         "classification": CLASSES[int(label)] if keep else "unresolved"})
    return rows


def calibrate(name, datasets, replicates, seed):
    chromosomes = geometry(name, seed)
    length = 250 if name.startswith("clustered") else 0
    truth = target(chromosomes, length)
    covered = np.zeros((len(BLOCKS), len(EDGES) - 1, 4), dtype=int)
    available = np.zeros_like(covered)
    widths = np.zeros_like(covered, dtype=float)
    estimates = []
    for iteration in range(datasets):
        rng = np.random.default_rng(np.random.SeedSequence([seed, iteration, 1]))
        rows = simulate(chromosomes, length, rng)
        result, _ = analyze(rows, EDGES)
        estimates.append(np.asarray([r["excess"] for r in result["reference_estimates"]]).reshape(truth.shape))
        intervals, _ = uncertainty(rows, EDGES, BLOCKS, replicates, seed + iteration)
        for interval in intervals:
            a = BLOCKS.index(interval["requested_block_bp"])
            b = EDGES.index(interval["distance_start"])
            k = METRICS.index(interval["metric"])
            if interval["status"] == "available":
                available[a, b, k] += 1
                widths[a, b, k] += interval["upper"] - interval["lower"]
                covered[a, b, k] += interval["lower"] <= truth[b, k] <= interval["upper"]
        if (iteration + 1) % 50 == 0:
            print(f"{name}: {iteration + 1}/{datasets}", flush=True)
    output = []
    bias = np.mean(estimates, axis=0) - truth
    for a, block in enumerate(BLOCKS):
        for b in range(len(EDGES) - 1):
            for k, metric in enumerate(METRICS):
                n, count = int(available[a, b, k]), int(covered[a, b, k])
                ci = binomtest(count, n).proportion_ci() if n else None
                coverage = count / n if n else None
                output.append({"block_bp": block, "distance_start": EDGES[b], "metric": metric,
                               "target": float(truth[b, k]), "bias": float(bias[b, k]),
                               "available": n, "withheld": datasets - n, "coverage": coverage,
                               "coverage_mc_interval": [ci.low, ci.high] if ci else None,
                               "mean_width": float(widths[a, b, k] / n) if n else None,
                               "passes_prespecified_coverage": bool(n >= 0.95 * datasets and
                                                                     0.90 <= count / n <= 0.99) if n else False})
    return {"scenario": name, "assumption_violation": name in ("frequency_gradient", "spatial_missingness"),
            "results": output}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datasets", type=int, default=200)
    parser.add_argument("--replicates", type=int, default=499)
    parser.add_argument("--seed", type=int, default=20260911)
    parser.add_argument("--scenarios", nargs="+", choices=SCENARIOS, default=list(SCENARIOS))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    result = {"datasets": args.datasets, "replicates": args.replicates, "seed": args.seed,
              "python": platform.python_version(), "numpy": np.__version__,
              "edges": EDGES, "blocks": BLOCKS, "scenarios": []}
    for name in args.scenarios:
        result["scenarios"].append(calibrate(name, args.datasets, args.replicates,
                                             args.seed + SCENARIOS.index(name) * 10000))
        result["seconds"] = time.perf_counter() - started
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    stationary = [r for s in result["scenarios"] if not s["assumption_violation"] for r in s["results"]]
    return 0 if all(r["passes_prespecified_coverage"] for r in stationary) else 1


if __name__ == "__main__":
    raise SystemExit(main())
