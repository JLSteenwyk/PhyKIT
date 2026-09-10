import json
from pathlib import Path

import pytest

from phykit.phykit import Phykit


FIXTURE = Path(__file__).resolve().parents[3] / "sample_files" / "topology_landscape"


def arguments(tmp_path, *extra):
    return ["--manifest", str(FIXTURE / "genes.tsv"),
            "--groups", str(FIXTURE / "groups.json"),
            "--coordinates", "reference", "bed", str(FIXTURE / "reference.bed"),
            "--output-prefix", str(tmp_path / "map"), "--window-bp", "500", "--json", *extra]


def test_end_to_end_counts_and_files(tmp_path, capsys):
    Phykit.topology_landscape(arguments(tmp_path, "--outgroup", "Outgroup"))
    result = json.loads(capsys.readouterr().out)
    assert [r["classification"] for r in result["genes"]] == [
        "topology_1", "topology_2", "topology_3", "unresolved",
        "insufficient_sampling", "incompatible_groups",
    ]
    assert result["overall"][0]["informative_genes"] == 3
    assert result["overall"][0]["total_genes"] == 6
    assert result["neighborhoods"][0]["topology_1_proportion"] == pytest.approx(1 / 3)
    assert result["neighborhoods"][1]["total_genes"] == 0
    assert result["labels"]["topology_1"] == "Ctenophores-sister"
    for output in result["output_files"].values():
        assert Path(output).is_file()


def test_reference_branch_and_gene_windows(tmp_path, capsys):
    args = arguments(tmp_path)
    index = args.index("--groups")
    del args[index:index + 2]
    index = args.index("--window-bp")
    del args[index:index + 2]
    args += ["--reference-tree", str(FIXTURE / "ctenophore.tre"), "--branch-taxa",
             "outgroup", "ctenophore", "--window-genes", "2"]
    Phykit.topology_landscape(args)
    result = json.loads(capsys.readouterr().out)
    assert [r["total_genes"] for r in result["neighborhoods"]] == [2, 2, 1, 1]
    assert list(result["groups"]) == ["A", "B", "C", "D"]


def test_multiple_references_and_plot(tmp_path, capsys, monkeypatch):
    monkeypatch.setenv("MPLBACKEND", "Agg")
    Phykit.topology_landscape(arguments(
        tmp_path, "--coordinates", "second", "bed", str(FIXTURE / "reference.bed"),
        "--plot", "--dpi", "60", "--labels", "Ctenophore-sister", "Sponge-sister", "Ctenophore + sponge",
    ))
    result = json.loads(capsys.readouterr().out)
    assert len(result["plots"]) == 4
    assert len(result["overall"]) == 2
    import matplotlib.image as mpimg
    import numpy as np
    for plot in result["plots"]:
        pixels = mpimg.imread(plot["path"])
        assert pixels.shape[0] > 100 and pixels.shape[1] > 100
        assert np.std(pixels[:, :, :3]) > 0.02


@pytest.mark.parametrize("options", [
    ["--min-support", "101"], ["--window-bp", "0"], ["--outgroup", "unknown"],
    ["--interval", "10", "5"], ["--chromosome", "absent"],
    ["--plot", "--plot-output", "bad.gif"], ["--plot", "--colors", "not_a_color"],
    ["--plot", "--fig-width", "nan"], ["--branch-taxa", "outgroup"],
])
def test_cli_errors(tmp_path, capsys, options):
    with pytest.raises(SystemExit) as error:
        Phykit.topology_landscape(arguments(tmp_path, *options))
    assert error.value.code == 2
    assert capsys.readouterr().out


def test_missing_coordinates_and_extra_ids(tmp_path, capsys):
    bed = tmp_path / "subset.bed"
    bed.write_text("chr1\t0\t10\tgene1\nchr1\t20\t30\tnot_a_tree\n")
    args = arguments(tmp_path)
    args[args.index(str(FIXTURE / "reference.bed"))] = str(bed)
    with pytest.raises(SystemExit):
        Phykit.topology_landscape(args)
    capsys.readouterr()
    Phykit.topology_landscape(args + ["--unmapped", "skip"])
    result = json.loads(capsys.readouterr().out)
    assert len(result["diagnostics"]["unmapped_genes"]) == 5
    assert result["diagnostics"]["coordinate_genes_without_trees"] == ["not_a_tree"]
    assert len(result["genes"]) == 1


@pytest.mark.parametrize("content", ["[]", "invalid", '{"A":["a"],"A":["b"]}'])
def test_bad_groups_file(tmp_path, capsys, content):
    groups = tmp_path / "groups.json"
    groups.write_text(content)
    args = arguments(tmp_path)
    args[args.index(str(FIXTURE / "groups.json"))] = str(groups)
    with pytest.raises(SystemExit) as error:
        Phykit.topology_landscape(args)
    assert error.value.code == 2
