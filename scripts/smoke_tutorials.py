#!/usr/bin/env python3
"""Run one portable command from every tutorial in an isolated directory."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import tarfile
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs" / "data"
MANIFEST = ROOT / "docs" / "_data" / "tutorial_smoke.json"


def prepare_data(destination: Path) -> None:
    for source in DATA.iterdir():
        if source.is_file():
            shutil.copy2(source, destination / source.name)
    for archive in destination.glob("*.tar.gz"):
        with tarfile.open(archive) as handle:
            root = destination.resolve()
            for member in handle.getmembers():
                target = (destination / member.name).resolve()
                if root not in (target, *target.parents):
                    raise RuntimeError(f"archive path escapes destination: {member.name}")
            handle.extractall(destination)


def run_smoke(python: str) -> int:
    tutorials = json.loads(MANIFEST.read_text())["tutorials"]
    failures = []
    with tempfile.TemporaryDirectory(prefix="phykit-docs-") as temporary:
        root = Path(temporary)
        prepare_data(root)
        environment = dict(os.environ, MPLBACKEND="Agg")
        for tutorial in tutorials:
            cwd = root / tutorial.get("cwd", "")
            command = [python, "-m", "phykit", *tutorial["command"]]
            result = subprocess.run(
                command,
                cwd=cwd,
                env=environment,
                text=True,
                capture_output=True,
                timeout=180,
            )
            missing = [
                artifact for artifact in tutorial.get("artifacts", [])
                if not (cwd / artifact).exists()
            ]
            if result.returncode or missing:
                failures.append(
                    {
                        "tutorial": tutorial["number"],
                        "command": command,
                        "returncode": result.returncode,
                        "missing": missing,
                        "stdout": result.stdout[-2000:],
                        "stderr": result.stderr[-2000:],
                    }
                )
            else:
                print(f"tutorial {tutorial['number']:02d}: passed")
    if failures:
        print(json.dumps(failures, indent=2))
        return 1
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--python", default=os.environ.get("PYTHON", "python"))
    args = parser.parse_args()
    interpreter = Path(args.python)
    if interpreter.exists():
        args.python = str(interpreter.resolve())
    return run_smoke(args.python)


if __name__ == "__main__":
    raise SystemExit(main())
