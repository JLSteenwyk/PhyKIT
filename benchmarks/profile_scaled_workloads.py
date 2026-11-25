#!/usr/bin/env python3
"""Time deterministic CLI inputs and optionally capture Python call profiles."""

import argparse
import json
from pathlib import Path
import random
import shutil
import statistics
import subprocess
import sys
import tempfile
import time


ROOT = Path(__file__).resolve().parents[1]
SIZES = {"small": (32, 300), "large": (1024, 1000)}
COMMANDS = ("pairwise_identity", "patristic_distances", "phylogenetic_signal")


def balanced_newick(names):
    if len(names) == 1:
        return f"{names[0]}:0.1"
    midpoint = len(names) // 2
    return f"({balanced_newick(names[:midpoint])},{balanced_newick(names[midpoint:])}):0.1"


def make_inputs(directory, taxa, sites):
    rng = random.Random(41)
    names = [f"taxon{index:04d}" for index in range(taxa)]
    alignment = directory / "alignment.fa"
    alignment.write_text(
        "".join(f">{name}\n{''.join(rng.choices('ACGT-', k=sites))}\n" for name in names)
    )
    tree = directory / "tree.nwk"
    tree.write_text(balanced_newick(names) + ";\n")
    traits = directory / "traits.tsv"
    traits.write_text("".join(f"{name}\t{rng.random():.12f}\n" for name in names))
    return {
        "pairwise_identity": [str(alignment)],
        "patristic_distances": [str(tree)],
        "phylogenetic_signal": ["-t", str(tree), "-d", str(traits), "-p", "100"],
    }


def measure(command, directory, repeat, gnu_time):
    subprocess.run(command, stdout=subprocess.DEVNULL, check=True)
    durations, peaks = [], []
    rss = directory / "rss.txt"
    for _ in range(repeat):
        measured = [gnu_time, "-f", "%M", "-o", str(rss), *command] if gnu_time else command
        started = time.perf_counter()
        subprocess.run(measured, stdout=subprocess.DEVNULL, check=True)
        durations.append(time.perf_counter() - started)
        if gnu_time:
            peaks.append(int(rss.read_text().strip()))
    return {
        "median_s": statistics.median(durations),
        "min_s": min(durations),
        "max_s": max(durations),
        "median_peak_rss_kib": statistics.median(peaks) if peaks else None,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeat", type=int, default=5)
    parser.add_argument("--size", choices=SIZES, action="append")
    parser.add_argument("--command", choices=COMMANDS, action="append")
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 4])
    parser.add_argument("--gnu-time", help="Path to GNU time, such as /usr/bin/time")
    parser.add_argument("--profile-dir", type=Path, help="Save a serial cProfile run per workload")
    args = parser.parse_args()
    if args.repeat < 1 or min(args.threads) < 1:
        parser.error("--repeat and --threads must be positive")
    if args.gnu_time and not shutil.which(args.gnu_time):
        parser.error("--gnu-time must be an executable")
    if args.profile_dir:
        args.profile_dir.mkdir(parents=True, exist_ok=True)

    for size in args.size or SIZES:
        taxa, sites = SIZES[size]
        with tempfile.TemporaryDirectory(prefix="phykit-scaled-") as scratch:
            directory = Path(scratch)
            workloads = make_inputs(directory, taxa, sites)
            for name in args.command or COMMANDS:
                for threads in args.threads:
                    cli = [str(ROOT / "phykit-runner.py"), "--threads", str(threads), name, *workloads[name]]
                    results = measure([sys.executable, *cli], directory, args.repeat, args.gnu_time)
                    print(json.dumps({
                        "command": name, "size": size, "taxa": taxa, "alignment_sites": sites,
                        "threads": threads, "repeat": args.repeat, **results,
                    }), flush=True)
                    if args.profile_dir and threads == 1:
                        profile = args.profile_dir / f"{name}-{size}.prof"
                        subprocess.run(
                            [sys.executable, "-m", "cProfile", "-o", str(profile), *cli],
                            stdout=subprocess.DEVNULL, check=True,
                        )


if __name__ == "__main__":
    main()
