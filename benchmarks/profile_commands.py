#!/usr/bin/env python3
"""Profile representative PhyKIT commands from a source checkout.

This script is intentionally small and dependency-free. It times real CLI
dispatch through ``phykit-runner.py`` against existing sample files, then sorts
commands by mean wall time so optimization passes can start from measured
bottlenecks instead of one-off ad hoc snippets.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import statistics
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "phykit-runner.py"
SAMPLES = ROOT / "tests" / "sample_files"

COMMANDS: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("alignment_length", ("alignment_length", str(SAMPLES / "test_alignment_0.fa"))),
    ("gc_content", ("gc_content", str(SAMPLES / "test_alignment_0.fa"))),
    (
        "gc_content_verbose",
        ("gc_content", str(SAMPLES / "test_alignment_0.fa"), "--verbose"),
    ),
    ("variable_sites", ("variable_sites", str(SAMPLES / "test_alignment_0.fa"))),
    (
        "composition_per_taxon",
        ("composition_per_taxon", str(SAMPLES / "test_alignment_0.fa")),
    ),
    ("alignment_entropy", ("alignment_entropy", str(SAMPLES / "test_alignment_0.fa"))),
    ("pairwise_identity", ("pairwise_identity", str(SAMPLES / "test_alignment_0.fa"))),
    (
        "pairwise_identity_verbose",
        ("pairwise_identity", str(SAMPLES / "test_alignment_0.fa"), "--verbose"),
    ),
    (
        "parsimony_informative_sites",
        ("parsimony_informative_sites", str(SAMPLES / "test_alignment_0.fa")),
    ),
    (
        "evolutionary_rate_per_site",
        ("evolutionary_rate_per_site", str(SAMPLES / "test_alignment_0.fa")),
    ),
    (
        "column_score",
        (
            "column_score",
            str(SAMPLES / "simple.fa"),
            "--reference",
            str(SAMPLES / "simple_reference.fa"),
        ),
    ),
    (
        "sum_of_pairs_score",
        (
            "sum_of_pairs_score",
            str(SAMPLES / "simple.fa"),
            "--reference",
            str(SAMPLES / "simple_reference.fa"),
        ),
    ),
    ("alignment_length_no_gaps", ("alignment_length_no_gaps", str(SAMPLES / "simple.fa"))),
    ("tip_labels", ("tip_labels", str(SAMPLES / "tree_simple.tre"))),
    ("dvmc", ("dvmc", str(SAMPLES / "tree_simple.tre"))),
    ("treeness", ("treeness", str(SAMPLES / "Yeasts_2832_eMRC_reference_renamed.tree"))),
    ("total_tree_length", ("total_tree_length", str(SAMPLES / "tree_simple.tre"))),
    (
        "evolutionary_rate",
        (
            "evolutionary_rate",
            str(SAMPLES / "12_YPR191W_Anc_7.548_codon_aln.fasta.clipkit.treefile"),
        ),
    ),
    (
        "bipartition_support_stats",
        ("bipartition_support_stats", str(SAMPLES / "small_Aspergillus_tree.tre")),
    ),
    (
        "internal_branch_stats",
        ("internal_branch_stats", str(SAMPLES / "small_Aspergillus_tree.tre")),
    ),
    (
        "monophyly_check_true",
        (
            "monophyly_check",
            str(SAMPLES / "small_Aspergillus_tree.tre"),
            str(SAMPLES / "small_Aspergillus_tree.monophyly_check.true.txt"),
        ),
    ),
    (
        "monophyly_check_false",
        (
            "monophyly_check",
            str(SAMPLES / "small_Aspergillus_tree.tre"),
            str(SAMPLES / "small_Aspergillus_tree.monophyly_check.false.txt"),
        ),
    ),
    ("patristic_distances", ("patristic_distances", str(SAMPLES / "tree_simple.tre"))),
    ("long_branch_score", ("long_branch_score", str(SAMPLES / "small_Aspergillus_tree.tre"))),
    (
        "long_branch_score_verbose",
        ("long_branch_score", str(SAMPLES / "small_Aspergillus_tree.tre"), "--verbose"),
    ),
    (
        "terminal_branch_stats",
        ("terminal_branch_stats", str(SAMPLES / "small_Aspergillus_tree.tre")),
    ),
    (
        "compositional_bias_per_site",
        ("compositional_bias_per_site", str(SAMPLES / "simple.fa")),
    ),
)


def _command_map() -> dict[str, tuple[str, ...]]:
    return dict(COMMANDS)


def _run_once(command_args: tuple[str, ...], env: dict[str, str]) -> float:
    start = time.perf_counter()
    subprocess.run(
        (sys.executable, str(RUNNER), *command_args),
        cwd=ROOT,
        env=env,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True,
    )
    return time.perf_counter() - start


def profile_commands(
    selected: list[str],
    repeat: int,
    warmup: int,
) -> list[dict[str, object]]:
    commands = _command_map()
    env = dict(os.environ)
    env["PYTHONPATH"] = str(ROOT)
    results: list[dict[str, object]] = []

    for name in selected:
        command_args = commands[name]
        for _ in range(warmup):
            _run_once(command_args, env)

        timings = [_run_once(command_args, env) for _ in range(repeat)]
        results.append(
            {
                "command": name,
                "mean_seconds": statistics.mean(timings),
                "median_seconds": statistics.median(timings),
                "min_seconds": min(timings),
                "runs": repeat,
            }
        )

    results.sort(key=lambda row: row["mean_seconds"], reverse=True)
    return results


def _parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Time representative PhyKIT commands and rank slowest means.",
    )
    parser.add_argument(
        "--repeat",
        type=int,
        default=3,
        help="Measured runs per command after warmup.",
    )
    parser.add_argument(
        "--warmup",
        type=int,
        default=1,
        help="Unmeasured warmup runs per command.",
    )
    parser.add_argument(
        "--command",
        action="append",
        choices=tuple(name for name, _ in COMMANDS),
        help="Command to profile. May be passed more than once.",
    )
    parser.add_argument("--json", action="store_true", help="Emit JSON rows.")
    args = parser.parse_args(argv)
    if args.repeat < 1:
        parser.error("--repeat must be at least 1")
    if args.warmup < 0:
        parser.error("--warmup must be at least 0")
    return args


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(sys.argv[1:] if argv is None else argv)
    selected = args.command or [name for name, _ in COMMANDS]
    results = profile_commands(selected, args.repeat, args.warmup)

    if args.json:
        print(json.dumps(results, indent=2))
        return 0

    print("command\tmean_s\tmedian_s\tmin_s\truns")
    for row in results:
        print(
            f"{row['command']}\t"
            f"{row['mean_seconds']:.6f}\t"
            f"{row['median_seconds']:.6f}\t"
            f"{row['min_seconds']:.6f}\t"
            f"{row['runs']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
