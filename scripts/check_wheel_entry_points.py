#!/usr/bin/env python3
"""Verify that a built wheel exposes every registered console script."""

from __future__ import annotations

import argparse
import configparser
import io
from pathlib import Path
from zipfile import ZipFile

from phykit.cli_registry import PUBLIC_COMMAND_TO_HANDLER


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("wheel", type=Path)
    args = parser.parse_args()
    with ZipFile(args.wheel) as archive:
        entry_points_path = next(
            name for name in archive.namelist() if name.endswith("entry_points.txt")
        )
        config = configparser.ConfigParser()
        config.optionxform = str
        config.read_file(io.StringIO(archive.read(entry_points_path).decode()))
    actual = dict(config["console_scripts"])
    expected = {"phykit": "phykit.phykit:main"}
    expected.update(
        {
            f"pk_{command}": f"phykit.phykit:{handler}"
            for command, handler in PUBLIC_COMMAND_TO_HANDLER.items()
        }
    )
    if actual != expected:
        missing = sorted(set(expected) - set(actual))
        unexpected = sorted(set(actual) - set(expected))
        wrong = sorted(
            name for name in set(actual) & set(expected)
            if actual[name] != expected[name]
        )
        print(f"missing entry points: {missing}")
        print(f"unexpected entry points: {unexpected}")
        print(f"incorrect handlers: {wrong}")
        return 1
    print(f"verified {len(actual)} console entry points in {args.wheel.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
