import json
import subprocess
import sys
from pathlib import Path

import pytest


@pytest.mark.integration
@pytest.mark.parametrize("command", ["topology_autocorrelation", "topo_ac"])
def test_autocorrelation_cli(tmp_path, command):
    fixture = Path(__file__).resolve().parents[1] / "sample_files" / "topology_landscape"
    result = subprocess.run([
        sys.executable, "-m", "phykit", command, "--manifest", str(fixture / "genes.tsv"),
        "--groups", str(fixture / "groups.json"), "--coordinates", "ref", "bed", str(fixture / "reference.bed"),
        "--output-prefix", str(tmp_path / "result"), "--distance-edges", "0", "101", "201", "--json",
    ], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    data = json.loads(result.stdout)
    assert [r["pair_count"] for r in data["reference_estimates"][::4]] == [2, 1]
    assert data["parameters"]["maximum_distance_exclusive"] == 201
