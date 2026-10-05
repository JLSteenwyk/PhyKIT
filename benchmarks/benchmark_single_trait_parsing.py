#!/usr/bin/env python3
"""Compare current service parsers with their bodies from a trusted Git revision."""

import argparse
import ast
from contextlib import redirect_stderr
import io
import json
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import tempfile
import time
from types import MethodType


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from phykit.errors import PhykitUserError
from phykit.services.tree.cont_map import ContMap
from phykit.services.tree.phenogram import Phenogram


def original_parser(revision, module, class_name):
    path = f"phykit/services/tree/{module}.py"
    source = subprocess.check_output(
        ["git", "show", f"{revision}:{path}"], cwd=ROOT, text=True,
    )
    tree = ast.parse(source)
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == class_name)
    method = next(n for n in cls.body if isinstance(n, ast.FunctionDef)
                  and n.name == "_parse_single_trait_data")
    namespace = {"sys": sys, "PhykitUserError": PhykitUserError}
    # Load only the original, state-independent parser, not the historical service.
    exec(compile(ast.Module(body=[method], type_ignores=[]), path, "exec"), namespace)
    return MethodType(namespace["_parse_single_trait_data"], object())


def measure(parser, path, tips, loops):
    with redirect_stderr(io.StringIO()):
        start = time.perf_counter()
        for _ in range(loops):
            parser(path, tips)
        return (time.perf_counter() - start) / loops


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-ref", required=True, help="trusted pre-extraction Git revision")
    parser.add_argument("--rows", type=int, default=500_000)
    parser.add_argument("--repeats", type=int, default=9)
    args = parser.parse_args()
    if args.rows < 4 or args.repeats < 1:
        parser.error("rows must be at least 4 and repeats must be positive")
    revision = subprocess.check_output(
        ["git", "rev-parse", "--verify", args.baseline_ref], cwd=ROOT, text=True,
    ).strip()
    results = []
    with tempfile.TemporaryDirectory(prefix="phykit-trait-parser-") as directory:
        large = Path(directory) / "large.tsv"
        names = [f"T{i}" for i in range(args.rows)]
        large.write_text("".join(f"{name}\t{i % 17}.5\n" for i, name in enumerate(names)))
        small = Path(directory) / "small.tsv"
        small.write_text("".join(f"T{i}\t{i}.5\n" for i in range(8)))
        cases = [
            ("small", str(small), names[:8] if args.rows >= 8 else [f"T{i}" for i in range(8)], 500),
            ("ordered", str(large), names, 1),
            ("reordered", str(large), list(reversed(names)), 1),
            ("partial_overlap", str(large), names[:-1] + ["missing"], 1),
        ]
        for module, cls in (("cont_map", ContMap), ("phenogram", Phenogram)):
            before = original_parser(revision, module, cls.__name__)
            after = cls.__new__(cls)._parse_single_trait_data
            for name, path, tips, loops in cases:
                old_warnings, new_warnings = io.StringIO(), io.StringIO()
                with redirect_stderr(old_warnings):
                    expected = before(path, tips)
                with redirect_stderr(new_warnings):
                    actual = after(path, tips)
                assert expected == actual and list(expected) == list(actual)
                assert old_warnings.getvalue() == new_warnings.getvalue()
                del expected, actual
                samples = {"before": [], "after": []}
                for repeat in range(args.repeats):
                    order = [("before", before), ("after", after)]
                    if repeat % 2:
                        order.reverse()
                    for label, function in order:
                        samples[label].append(measure(function, path, tips, loops))
                baseline = statistics.median(samples["before"])
                current = statistics.median(samples["after"])
                results.append({"service": module, "case": name,
                                "baseline_seconds": baseline, "current_seconds": current,
                                "baseline_over_current": baseline / current,
                                "samples": samples})
    print(json.dumps({"baseline_ref": revision, "python": platform.python_version(),
                      "platform": platform.platform(), "rows": args.rows,
                      "repeats": args.repeats, "results": results}, indent=2))


if __name__ == "__main__":
    main()
