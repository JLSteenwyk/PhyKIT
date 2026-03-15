"""
Unit tests for phylo_heatmap (phylogenetic heatmap visualization).

Tests initialization, trait matrix parsing, plot creation, and
standardization. Analogous to R's phytools::phylo.heatmap().
"""
import pytest
from argparse import Namespace
from pathlib import Path

from phykit.services.tree.phylo_heatmap import PhyloHeatmap


@pytest.fixture
def args(tmp_path):
    return Namespace(
        tree="tests/sample_files/tree_simple.tre",
        data="tests/sample_files/tree_simple_multi_traits.tsv",
        output=str(tmp_path / "heatmap.png"),
        split=0.3,
        standardize=False,
        cmap="viridis",
    )


class TestPhyloHeatmapInit:
    def test_init_sets_fields(self, args):
        ph = PhyloHeatmap(args)
        assert ph.tree_file_path == args.tree
        assert ph.data_path == args.data
        assert ph.split == 0.3
        assert ph.standardize is False
        assert ph.cmap_name == "viridis"
        assert ph.json_output is False

    def test_invalid_split_exits(self, tmp_path):
        args = Namespace(
            tree="t.tre", data="d.tsv",
            output=str(tmp_path / "out.png"),
            split=0.0,
        )
        with pytest.raises(SystemExit):
            PhyloHeatmap(args)

    def test_invalid_split_above_one_exits(self, tmp_path):
        args = Namespace(
            tree="t.tre", data="d.tsv",
            output=str(tmp_path / "out.png"),
            split=1.0,
        )
        with pytest.raises(SystemExit):
            PhyloHeatmap(args)


class TestTraitParsing:
    def test_parse_multi_column_tsv(self, args):
        ph = PhyloHeatmap(args)
        tree = ph.read_tree_file()
        tree_tips = [t.name for t in tree.get_terminals()]
        trait_names, trait_data = ph._parse_trait_matrix(args.data, tree_tips)
        assert trait_names == ["body_mass", "brain_size", "longevity"]
        assert len(trait_data) == 8
        assert "raccoon" in trait_data
        assert len(trait_data["raccoon"]) == 3

    def test_missing_file_exits(self, args):
        ph = PhyloHeatmap(args)
        with pytest.raises(SystemExit):
            ph._parse_trait_matrix("/nonexistent.tsv", ["A", "B", "C"])

    def test_too_few_shared_taxa_exits(self, args, tmp_path):
        ph = PhyloHeatmap(args)
        f = tmp_path / "traits.tsv"
        f.write_text("taxon\tval\nX\t1.0\nY\t2.0\n")
        with pytest.raises(SystemExit):
            ph._parse_trait_matrix(str(f), ["A", "B", "C"])

    def test_non_numeric_exits(self, args, tmp_path):
        ph = PhyloHeatmap(args)
        f = tmp_path / "traits.tsv"
        f.write_text("taxon\tval\nA\tabc\nB\t2.0\nC\t3.0\n")
        with pytest.raises(SystemExit):
            ph._parse_trait_matrix(str(f), ["A", "B", "C"])


class TestPhyloHeatmapPlot:
    def test_creates_png(self, args):
        ph = PhyloHeatmap(args)
        ph.run()
        assert Path(args.output).exists()
        assert Path(args.output).stat().st_size > 0

    def test_creates_pdf(self, tmp_path):
        args = Namespace(
            tree="tests/sample_files/tree_simple.tre",
            data="tests/sample_files/tree_simple_multi_traits.tsv",
            output=str(tmp_path / "heatmap.pdf"),
            split=0.3, standardize=False, cmap="viridis",
        )
        ph = PhyloHeatmap(args)
        ph.run()
        assert Path(args.output).exists()

    def test_standardize_flag(self, tmp_path):
        args = Namespace(
            tree="tests/sample_files/tree_simple.tre",
            data="tests/sample_files/tree_simple_multi_traits.tsv",
            output=str(tmp_path / "heatmap_std.png"),
            split=0.3, standardize=True, cmap="viridis",
        )
        ph = PhyloHeatmap(args)
        ph.run()
        assert Path(args.output).exists()

    def test_custom_cmap(self, tmp_path):
        args = Namespace(
            tree="tests/sample_files/tree_simple.tre",
            data="tests/sample_files/tree_simple_multi_traits.tsv",
            output=str(tmp_path / "heatmap_hot.png"),
            split=0.3, standardize=False, cmap="hot",
        )
        ph = PhyloHeatmap(args)
        ph.run()
        assert Path(args.output).exists()

    def test_custom_split(self, tmp_path):
        args = Namespace(
            tree="tests/sample_files/tree_simple.tre",
            data="tests/sample_files/tree_simple_multi_traits.tsv",
            output=str(tmp_path / "heatmap_split.png"),
            split=0.5, standardize=False, cmap="viridis",
        )
        ph = PhyloHeatmap(args)
        ph.run()
        assert Path(args.output).exists()


class TestPhyloHeatmapJson:
    def test_json_output(self, mocker, args):
        args.json = True
        ph = PhyloHeatmap(args)
        mocked_json = mocker.patch("phykit.services.tree.phylo_heatmap.print_json")
        ph.run()
        payload = mocked_json.call_args.args[0]
        assert payload["n_taxa"] == 8
        assert payload["n_traits"] == 3
        assert payload["trait_names"] == ["body_mass", "brain_size", "longevity"]
        assert len(payload["tip_order"]) == 8
