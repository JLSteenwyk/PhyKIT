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
import sys
import tempfile
from zipfile import ZipFile


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from phykit.cli_registry import (  # noqa: E402
    COMMAND_IDENTITIES,
    PUBLIC_COMMAND_TO_HANDLER,
)


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
    python = python.absolute()
    with tempfile.TemporaryDirectory() as working_directory:
        completed = subprocess.run(
            [str(python), "-c", INSTALLED_TARGET_CHECK],
            check=False,
            capture_output=True,
            cwd=working_directory,
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


def resolve_wheel(path: Path) -> Path:
    """Resolve a wheel path or a directory containing exactly one wheel."""
    if path.is_file() and path.suffix == ".whl":
        return path
    if path.is_dir():
        wheels = sorted(path.glob("*.whl"))
        if len(wheels) == 1:
            return wheels[0]
        raise ValueError(f"expected one wheel in {path}, found {len(wheels)}")
    raise ValueError(f"wheel does not exist: {path}")


def create_installed_wheel_environment(wheel: Path, environment: Path) -> Path:
    """Install a wheel into a new environment and return its interpreter."""
    subprocess.run(
        [sys.executable, "-m", "venv", str(environment)],
        check=True,
    )
    if sys.platform == "win32":
        python = environment / "Scripts" / "python.exe"
    else:
        python = environment / "bin" / "python"
    subprocess.run(
        [
            str(python),
            "-m",
            "pip",
            "install",
            "--disable-pip-version-check",
            str(wheel.resolve()),
        ],
        check=True,
    )
    return python


def smoke_installed_targets(python: Path) -> dict[str, str]:
    """Return all metadata-load and canonical-help smoke failures."""
    errors = {
        f"installed target {name}": error
        for name, error in installed_target_load_errors(python).items()
    }
    errors.update(
        {
            f"installed help {name}": error
            for name, error in canonical_help_errors(python).items()
        }
    )
    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("wheel", type=Path, help="wheel file or containing directory")
    smoke_group = parser.add_mutually_exclusive_group()
    smoke_group.add_argument(
        "--smoke-python",
        type=Path,
        help="interpreter from an environment where the wheel is installed",
    )
    smoke_group.add_argument(
        "--install-smoke",
        action="store_true",
        help="install the wheel in a temporary environment before smoke testing",
    )
    args = parser.parse_args()
    try:
        wheel = resolve_wheel(args.wheel)
    except ValueError as error:
        parser.error(str(error))
    with ZipFile(wheel) as archive:
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
    if not args.smoke_python and not args.install_smoke:
        load_errors = target_load_errors(actual)
        if load_errors:
            for name, error in sorted(load_errors.items()):
                print(f"unloadable entry point {name}: {error}")
            return 1
    smoke_python = args.smoke_python
    temporary_environment = None
    if args.install_smoke:
        temporary_environment = tempfile.TemporaryDirectory()
        try:
            smoke_python = create_installed_wheel_environment(
                wheel, Path(temporary_environment.name) / "wheel-smoke"
            )
        except subprocess.CalledProcessError as error:
            temporary_environment.cleanup()
            print(f"failed to install wheel for smoke test: {error}")
            return 1
    if smoke_python:
        smoke_errors = smoke_installed_targets(smoke_python)
        if temporary_environment is not None:
            temporary_environment.cleanup()
        if smoke_errors:
            for name, error in sorted(smoke_errors.items()):
                print(f"failed {name}: {error}")
            return 1
    print(
        f"verified {len(actual)} loadable console entry points "
        f"in {wheel.name}"
    )
    if smoke_python:
        print(
            f"verified installed targets and {len(COMMAND_IDENTITIES)} "
            "canonical help commands"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
