#!/usr/bin/env python3
"""Report which PhyKIT definitions were observed in cProfile captures."""

import argparse
import ast
import csv
from pathlib import Path
import pstats
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("profiles", nargs="+", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    observations = {}
    for profile in args.profiles:
        for filename, line, name in pstats.Stats(str(profile)).stats:
            source = Path(filename).resolve()
            if source.is_relative_to(root / "phykit"):
                key = (source, line, name)
                observations.setdefault(key, []).append(profile.name)

    writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    writer.writerow(("file", "line", "function", "profiles"))
    for path in sorted((root / "phykit").rglob("*.py")):
        nodes = ast.walk(ast.parse(path.read_bytes(), filename=str(path)))
        definitions = sorted(
            (node for node in nodes if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))),
            key=lambda node: node.lineno,
        )
        for node in definitions:
            compiled_line = min([node.lineno, *(item.lineno for item in node.decorator_list)])
            measured = observations.get((path, compiled_line, node.name), [])
            writer.writerow((path.relative_to(root), node.lineno, node.name, ",".join(sorted(measured))))


if __name__ == "__main__":
    main()
