import argparse
from statistics import median
from time import perf_counter

import numpy as np
from scipy.linalg import cho_factor, solve_triangular

from phykit.services.tree.phylogenetic_signal import (
    PhylogeneticSignal,
    _permutation_p_value_ge,
)


def solve_k(x, covariance, permutations):
    n = len(x)
    factor, lower = cho_factor(covariance, lower=True, check_finite=False)
    ones = solve_triangular(factor, np.ones(n), lower=lower, check_finite=False)
    weight = float(ones @ ones)
    expected = (np.trace(covariance) - n / weight) / (n - 1)

    def evaluate(values):
        whitened = solve_triangular(
            factor, values.T, lower=lower, check_finite=False
        ).T
        means = (whitened @ ones) / weight
        residual = values - means[:, None]
        whitened -= means[:, None] * ones
        return (
            np.einsum("ij,ij->i", residual, residual)
            / np.einsum("ij,ij->i", whitened, whitened)
            / expected
        )

    observed = float(evaluate(x[None, :])[0])
    rng = np.random.default_rng(seed=42)
    sampled = np.empty(permutations)
    for start in range(0, permutations, 128):
        stop = min(start + 128, permutations)
        samples = np.array([rng.permutation(x) for _ in range(start, stop)])
        sampled[start:stop] = evaluate(samples)
    return dict(K=observed, p_value=_permutation_p_value_ge(sampled, observed))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeat", type=int, default=5)
    args = parser.parse_args()
    current = object.__new__(PhylogeneticSignal)._blombergs_k
    rng = np.random.default_rng(41)
    for n in (128, 1024):
        # Shared ancestry in a balanced tree, with equal terminal depths.
        levels = np.bitwise_xor.outer(np.arange(n), np.arange(n))
        depth = np.log2(np.maximum(levels, 1)).astype(int)
        c = (n.bit_length() - 2 - depth) * 0.1
        np.fill_diagonal(c, (n.bit_length() - 1) * 0.1)
        x = rng.normal(size=n)
        for p in (100, 1000):
            vals = []
            for fn in (current, solve_k):
                result = fn(x, c, p)
                times = []
                for _ in range(args.repeat):
                    start = perf_counter()
                    fn(x, c, p)
                    times.append(perf_counter() - start)
                vals.append((round(median(times), 6), round(min(times), 6), round(max(times), 6), result))
            print(n, p, vals, flush=True)

    # Symmetric covariance and repeated values put permutations on ties.
    for name, c, x in [
        ("identity", np.eye(8), np.arange(8.0)),
        ("repeated", np.eye(8) + 0.01, np.array([0.0] * 4 + [1.0] * 4)),
        ("near-singular", np.full((8,8), 1.0) + 1e-9 * np.eye(8), np.arange(8.0)),
    ]:
        print(name, current(x, c, 1000), solve_k(x, c, 1000), flush=True)

if __name__ == "__main__":
    main()
