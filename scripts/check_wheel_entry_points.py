#!/usr/bin/env python3
"""Verify that a built wheel exposes every registered console script."""

from __future__ import annotations

import argparse
import configparser
import importlib
import io
from pathlib import Path
from zipfile import ZipFile

from phykit.cli_registry import PUBLIC_COMMAND_TO_HANDLER


def target_load_errors(entry_points: dict[str, str]) -> dict[str, str]:
    """Return import or callability errors for console-script targets."""
    errors = {}
    for name, target in entry_points.items():
        module_name, separator, attribute = target.partition(":")
        if not separator or not module_name or not attribute:
            errors[name] = f"invalid target: {target}"
            continue
        try:
            value = getattr(importlib.import_module(module_name), attribute)
        except (AttributeError, ImportError) as error:
            errors[name] = f"{target}: {error}"
            continue
        if not callable(value):
            errors[name] = f"target is not callable: {target}"
    return errors


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
    load_errors = target_load_errors(actual)
    if load_errors:
        for name, error in sorted(load_errors.items()):
            print(f"unloadable entry point {name}: {error}")
        return 1
    print(
        f"verified {len(actual)} loadable console entry points "
        f"in {args.wheel.name}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
