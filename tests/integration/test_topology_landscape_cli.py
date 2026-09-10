import json
import os
from pathlib import Path
import subprocess
import sys

import pytest


@pytest.mark.integration
@pytest.mark.parametrize("command", ["topology_landscape", "topomap"])
def test_cli_genomic_landscape(tmp_path, command):
    fixture = Path(__file__).resolve().parents[1] / "sample_files" / "topology_landscape"
    result = subprocess.run([
        sys.executable, "-m", "phykit", command,
        "--manifest", str(fixture / "genes.tsv"),
        "--groups", str(fixture / "groups.json"),
        "--coordinates", "ref", "bed", str(fixture / "reference.bed"),
        "--output-prefix", str(tmp_path / "map"), "--json",
    ], capture_output=True, text=True, env=dict(os.environ, MPLBACKEND="Agg"))
    assert result.returncode == 0, result.stdout + result.stderr
    payload = json.loads(result.stdout)
    assert payload["overall"][0]["total_genes"] == 6
    assert payload["overall"][0]["informative_genes"] == 3
