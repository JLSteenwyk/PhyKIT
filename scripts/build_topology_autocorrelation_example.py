"""Build the downloadable synthetic distance-decay tutorial data."""

import csv
import shutil
import tarfile
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def main():
    rng = np.random.default_rng(712)
    with tempfile.TemporaryDirectory() as temporary:
        directory = Path(temporary) / "topology_autocorrelation"
        shutil.copytree(ROOT / "tests/sample_files/topology_autocorrelation", directory)
        with open(directory / "simulated.tsv", "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["reference", "chromosome", "gene_id", "start", "end", "anchor", "classification"])
            for c, n in enumerate((6000, 4000), 1):
                positions = np.cumsum(rng.integers(25, 176, n))
                probabilities = [0.5, 0.3, 0.2] if c == 1 else [0.2, 0.3, 0.5]
                label, previous = rng.choice(3, p=probabilities), 0
                for i, x in enumerate(positions):
                    if rng.random() >= np.exp(-(x - previous) / 250):
                        label = rng.choice(3, p=probabilities)
                    category = "unresolved" if rng.random() < 0.2 else f"topology_{label + 1}"
                    writer.writerow(["simulated", f"chr{c}", f"c{c}g{i}", x, x + 1, x, category])
                    previous = x
        target = ROOT / "docs/data/topology_autocorrelation_tutorial.tar.gz"
        with tarfile.open(target, "w:gz") as archive:
            archive.add(directory, arcname=directory.name)
        print(target)


if __name__ == "__main__":
    main()
