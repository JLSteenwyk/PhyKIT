"""Measure time and peak Python/NumPy allocation without storing all gene pairs."""

import argparse
import json
import platform
import time
import tracemalloc
from pathlib import Path

import numpy as np

from phykit.helpers.topology_autocorrelation import analyze


def benchmark(n, maximum):
    rows = [{"reference": "ref", "chromosome": "chr1", "gene_id": f"g{i}",
             "start": i * 100, "end": i * 100 + 1, "anchor": i * 100,
             "classification": f"topology_{i % 3 + 1}"} for i in range(n)]
    edges = np.linspace(0, maximum, 11, dtype=int).tolist()
    tracemalloc.start()
    started = time.perf_counter()
    result, _ = analyze(rows, edges)
    elapsed = time.perf_counter() - started
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    pairs = sum(r["pair_count"] for r in result["reference_estimates"] if r["metric"] == "agreement")
    neighbors = min(n - 1, (maximum - 1) // 100)
    expected = neighbors * n - neighbors * (neighbors + 1) // 2
    assert pairs == expected
    return {"genes": n, "bins": 10, "maximum_distance_exclusive": maximum, "pairs": pairs,
            "seconds_with_tracemalloc": elapsed, "peak_additional_traced_mib": peak / 2**20}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", nargs="+", type=int, default=[1000, 10000, 100000])
    parser.add_argument("--maxima", nargs="+", type=int, default=[1000, 1000000])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data = {"python": platform.python_version(), "platform": platform.platform(), "numpy": np.__version__,
            "memory_scope": "tracemalloc peak after input construction, including validation copies and NumPy buffers",
            "results": []}
    for n in args.sizes:
        for maximum in args.maxima:
            result = benchmark(n, maximum)
            data["results"].append(result)
            print(json.dumps(result), flush=True)
    args.output.write_text(json.dumps(data, indent=2) + "\n")


if __name__ == "__main__":
    main()
