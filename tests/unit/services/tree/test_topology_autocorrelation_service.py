import csv
import json
from pathlib import Path

import numpy as np
import pytest

from phykit.phykit import Phykit

FIXTURE = Path(__file__).resolve().parents[3] / "sample_files" / "topology_landscape"


def arguments(tmp_path, *extra):
    return ["--manifest", str(FIXTURE / "genes.tsv"), "--groups", str(FIXTURE / "groups.json"),
            "--coordinates", "ref", "bed", str(FIXTURE / "reference.bed"),
            "--output-prefix", str(tmp_path / "map"), "--json", *extra]


def test_raw_and_classified_paths_match(tmp_path, capsys):
    Phykit.topology_autocorrelation(arguments(tmp_path, "--distance-edges", "0", "101", "201"))
    raw = json.loads(capsys.readouterr().out)
    assert [r["pair_count"] for r in raw["reference_estimates"][::4]] == [2, 1]
    assert all(r["excess"] == 0 for r in raw["reference_estimates"])
    assert raw["clustering_range"]["estimate"] is None
    assert raw["uncertainty"] == []
    Phykit.topology_autocorrelation(["--classified-genes", raw["output_files"]["genes"],
                                    "--output-prefix", str(tmp_path / "second"),
                                    "--distance-edges", "0", "101", "201", "--json"])
    classified = json.loads(capsys.readouterr().out)
    assert raw["reference_estimates"] == classified["reference_estimates"]
    assert raw["chromosome_estimates"] == classified["chromosome_estimates"]
    for path in classified["output_files"].values():
        assert Path(path).stat().st_size > 0
    on_disk = json.loads(Path(classified["output_files"]["json"]).read_text())
    assert on_disk == classified


def test_raw_reference_edge_and_selection(tmp_path, capsys):
    args = arguments(tmp_path, "--bin-width", "101", "--max-distance", "250",
                     "--chromosome", "chr1", "--interval", "100", "300")
    i = args.index("--groups")
    del args[i:i + 2]
    args += ["--reference-tree", str(FIXTURE / "ctenophore.tre"), "--branch-taxa", "outgroup", "ctenophore"]
    Phykit.topology_autocorrelation(args)
    result = json.loads(capsys.readouterr().out)
    assert len(result["genes"]) == 2
    assert result["parameters"]["distance_edges"] == [0, 101, 202, 250]
    assert result["input_metadata"]["groups"]


@pytest.mark.parametrize("extension", ["png", "pdf", "svg"])
def test_multireference_plots_and_withheld_intervals(tmp_path, capsys, extension, monkeypatch):
    monkeypatch.setenv("MPLBACKEND", "Agg")
    Phykit.topology_autocorrelation(arguments(
        tmp_path, "--coordinates", "second", "bed", str(FIXTURE / "reference.bed"),
        "--plot", "--plot-output", str(tmp_path / f"figure.{extension}"), "--dpi", "60",
        "--distance-edges", "0", "101", "201", "--block-sizes", "50", "100", "--replicates", "499",
        "--labels", "Ctenophore-sister", "Sponge-sister", "Ctenophore + sponge"))
    result = json.loads(capsys.readouterr().out)
    assert len(result["plots"]) == 2
    assert all(r["status"] == "withheld" for r in result["uncertainty"])
    assert result["block_diagnostics"]
    for path in result["plots"]:
        assert Path(path).stat().st_size > 1000
        if extension == "png":
            from matplotlib.image import imread
            pixels = imread(path)
            assert pixels.shape[0] > 100
            assert np.std(pixels[:, :, :3]) > 0.02


@pytest.mark.parametrize("options", [
    ["--distance-edges", "0", "10", "--max-distance", "30"],
    ["--distance-edges", "1", "10"], ["--bin-width", "0"], ["--max-distance", "0"],
    ["--bin-width", "1", "--max-distance", "2000"], ["--seed", "-1"],
    ["--block-sizes", "10"], ["--block-sizes", "10", "20", "--replicates", "10"],
    ["--interval", "10", "1"], ["--chromosome", "absent"],
    ["--plot", "--fig-width", "-1"], ["--plot", "--fig-height", "nan"],
    ["--plot", "--fig-width", "0"], ["--plot", "--fig-height", "0"],
    ["--plot", "--colors", "not_a_color"], ["--plot", "--plot-output", "x.bad"],
])
def test_user_errors(tmp_path, capsys, options):
    with pytest.raises(SystemExit):
        Phykit.topology_autocorrelation(arguments(tmp_path, *options))
    assert "Traceback" not in capsys.readouterr().err


def test_output_protection_and_raw_options_rejected_for_classified(tmp_path, capsys):
    Phykit.topology_autocorrelation(arguments(tmp_path))
    capsys.readouterr()
    original = (tmp_path / "map.genes.tsv").read_bytes()
    base = ["--classified-genes", str(tmp_path / "map.genes.tsv"),
            "--output-prefix", str(tmp_path / "map")]
    with pytest.raises(SystemExit):
        Phykit.topology_autocorrelation(base)
    assert "overwrite" in capsys.readouterr().out
    assert (tmp_path / "map.genes.tsv").read_bytes() == original
    with pytest.raises(SystemExit):
        Phykit.topology_autocorrelation(base + ["--min-support", "50"])
    assert "cannot be combined" in capsys.readouterr().out
    with pytest.raises(SystemExit):
        Phykit.topology_autocorrelation(arguments(tmp_path, "--plot", "--plot-output", str(tmp_path / "map.json")))
    with pytest.raises(SystemExit):
        Phykit.topology_autocorrelation(arguments(tmp_path, "--output-prefix", str(tmp_path / "missing" / "out")))


def test_malformed_classified_and_missing_raw_coordinates(tmp_path, capsys):
    path = tmp_path / "bad.tsv"
    path.write_text("gene_id\nfoo\n")
    with pytest.raises(SystemExit):
        Phykit.topology_autocorrelation(["--classified-genes", str(path), "-o", str(tmp_path / "out")])
    args = arguments(tmp_path)
    i = args.index("--coordinates")
    del args[i:i + 4]
    with pytest.raises(SystemExit):
        Phykit.topology_autocorrelation(args)


def test_support_filtering_and_plain_output(tmp_path, capsys):
    args = arguments(tmp_path, "--min-support", "100", "--missing-support", "collapse")
    args.remove("--json")
    Phykit.topology_autocorrelation(args)
    captured = capsys.readouterr()
    assert "excess_agreement" in captured.out
    with open(tmp_path / "map.autocorrelation.tsv") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert all(int(r["pair_count"]) == 0 for r in rows)
