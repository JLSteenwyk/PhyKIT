#!/usr/bin/env python3
"""Verify that a built wheel exposes every registered console script."""

from __future__ import annotations

import argparse
import configparser
import importlib
import io
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
from zipfile import ZipFile

from phykit.cli_registry import COMMAND_IDENTITIES, PUBLIC_COMMAND_TO_HANDLER


INSTALLED_TARGET_CHECK = """
import json
from importlib.metadata import distribution

errors = {}
entry_points = distribution("phykit").entry_points
for entry_point in entry_points:
    if entry_point.group != "console_scripts":
        continue
    try:
        target = entry_point.load()
    except Exception as error:
        errors[entry_point.name] = f"{type(error).__name__}: {error}"
        continue
    if not callable(target):
        errors[entry_point.name] = f"target is not callable: {entry_point.value}"
print(json.dumps(errors, sort_keys=True))
"""


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


def installed_target_load_errors(python: Path) -> dict[str, str]:
    """Load console-script targets from the interpreter's installed wheel."""
    completed = subprocess.run(
        [str(python), "-c", INSTALLED_TARGET_CHECK],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode:
        return {
            "<metadata>": completed.stderr.strip() or completed.stdout.strip()
        }
    try:
        return json.loads(completed.stdout)
    except json.JSONDecodeError as error:
        return {"<metadata>": f"invalid checker output: {error}"}


def canonical_help_errors(python: Path) -> dict[str, str]:
    """Run help through every installed canonical standalone executable."""
    scripts_directory = python.absolute().parent
    errors = {}
    with tempfile.TemporaryDirectory() as working_directory:
        for identity in COMMAND_IDENTITIES:
            name = f"pk_{identity.canonical}"
            executable = shutil.which(name, path=str(scripts_directory))
            if executable is None:
                errors[name] = f"not found in {scripts_directory}"
                continue
            completed = subprocess.run(
                [executable, "--help"],
                check=False,
                capture_output=True,
                cwd=working_directory,
                text=True,
            )
            if completed.returncode:
                detail = completed.stderr.strip() or completed.stdout.strip()
                errors[name] = f"exit {completed.returncode}: {detail}"
    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("wheel", type=Path)
    parser.add_argument(
        "--smoke-python",
        type=Path,
        help="interpreter from an environment where the wheel is installed",
    )
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
    if args.smoke_python:
        installed_errors = installed_target_load_errors(args.smoke_python)
        if installed_errors:
            for name, error in sorted(installed_errors.items()):
                print(f"unloadable installed entry point {name}: {error}")
            return 1
        help_errors = canonical_help_errors(args.smoke_python)
        if help_errors:
            for name, error in sorted(help_errors.items()):
                print(f"failed installed help smoke {name}: {error}")
            return 1
    print(
        f"verified {len(actual)} loadable console entry points "
        f"in {args.wheel.name}"
    )
    if args.smoke_python:
        print(
            f"verified installed targets and {len(COMMAND_IDENTITIES)} "
            "canonical help commands"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
