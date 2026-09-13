"""Installed-package checks, copied into the Bioconda recipe by sync.py."""

import importlib
from importlib.metadata import distribution
from pathlib import Path
import shutil
import subprocess
import tempfile


entries = [e for e in distribution("phykit").entry_points if e.group == "console_scripts"]
assert entries, "Missing console script metadata"
for entry in entries:
    executable = shutil.which(entry.name)
    assert executable, f"Missing executable: {entry.name}"
    assert callable(entry.load()), f"Invalid target: {entry.value}"
    result = subprocess.run(
        [executable, "--help"], capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 0, f"{entry.name}: {result.stdout}\n{result.stderr}"

for module in ("Bio", "matplotlib", "numpy", "scipy", "sklearn", "tqdm", "umap"):
    importlib.import_module(module)

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt

with tempfile.TemporaryDirectory() as directory:
    tree = Path(directory) / "tree.nwk"
    tree.write_text("((A:1,B:1):1,C:1);\n")
    result = subprocess.run(
        ["pk_total_tree_length", str(tree)],
        capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 0, result.stderr
    assert float(result.stdout.strip()) == 4.0, result.stdout
    figure, axes = plt.subplots()
    axes.plot([0, 1], [0, 1])
    output = Path(directory) / "plot.png"
    figure.savefig(output)
    plt.close(figure)
    assert output.stat().st_size > 0

print(f"Verified {len(entries)} console scripts, runtime imports, tree length, and plotting")
