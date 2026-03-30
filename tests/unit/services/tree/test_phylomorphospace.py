import os

import pytest
import json
import numpy as np
from argparse import Namespace
from pathlib import Path

from phykit.services.tree.phylomorphospace import Phylomorphospace
from phykit.helpers.trait_parsing import parse_multi_trait_file
from phykit.errors import PhykitUserError


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
MULTI_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_multi_traits.tsv")


@pytest.fixture
def default_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=MULTI_TRAITS_FILE,
        trait_x="body_mass",
        trait_y="brain_size",
        color_by=None,
        plot_output="phylomorphospace_plot.png",
        json=False,
    )


@pytest.fixture
def auto_select_two_trait_file(tmp_path):
    """Create a trait file with exactly 2 traits for auto-selection."""
    trait_file = tmp_path / "two_traits.tsv"
    trait_file.write_text(
        "taxon\ttrait_a\ttrait_b\n"
        "raccoon\t1.0\t2.0\n"
        "bear\t2.0\t3.0\n"
        "sea_lion\t3.0\t4.0\n"
        "seal\t4.0\t5.0\n"
        "monkey\t5.0\t6.0\n"
        "cat\t6.0\t7.0\n"
        "weasel\t7.0\t8.0\n"
        "dog\t8.0\t9.0\n"
    )
    return str(trait_file)


class TestProcessArgs:
    def test_defaults(self):
        args = Namespace(tree="t.tre", trait_data="d.tsv")
        svc = Phylomorphospace(args)
        assert svc.tree_file_path == "t.tre"
        assert svc.trait_data_path == "d.tsv"
        assert svc.trait_x is None
        assert svc.trait_y is None
        assert svc.color_by is None
        assert svc.plot_output == "phylomorphospace_plot.png"
        assert svc.json_output is False

    def test_overrides(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            trait_x="body_mass",
            trait_y="brain_size",
            color_by="longevity",
            plot_output="my_plot.png",
            json=True,
        )
        svc = Phylomorphospace(args)
        assert svc.trait_x == "body_mass"
        assert svc.trait_y == "brain_size"
        assert svc.color_by == "longevity"
        assert svc.plot_output == "my_plot.png"
        assert svc.json_output is True

    def test_auto_select_two_trait_file(self, auto_select_two_trait_file, tmp_path):
        plot_path = str(tmp_path / "auto.png")
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=auto_select_two_trait_file,
            trait_x=None,
            trait_y=None,
            color_by=None,
            plot_output=plot_path,
            json=False,
        )
        svc = Phylomorphospace(args)
        svc.run()
        assert Path(plot_path).exists()


class TestRun:
    def test_basic_text_output(self, default_args, tmp_path, capsys):
        default_args.plot_output = str(tmp_path / "test.png")
        svc = Phylomorphospace(default_args)
        svc.run()
        captured = capsys.readouterr()
        assert "Saved phylomorphospace plot" in captured.out
        assert Path(default_args.plot_output).exists()

    def test_json_output(self, default_args, tmp_path, capsys):
        default_args.json = True
        default_args.plot_output = str(tmp_path / "test.png")
        svc = Phylomorphospace(default_args)
        svc.run()
        captured = capsys.readouterr()
        payload = json.loads(captured.out)
        assert payload["trait_x"] == "body_mass"
        assert payload["trait_y"] == "brain_size"
        assert "tip_data" in payload
        assert "raccoon" in payload["tip_data"]
        assert "plot_output" in payload

    def test_trait_not_found_error(self, tmp_path):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            trait_x="nonexistent",
            trait_y="brain_size",
            color_by=None,
            plot_output=str(tmp_path / "test.png"),
            json=False,
        )
        svc = Phylomorphospace(args)
        with pytest.raises(PhykitUserError):
            svc.run()

    def test_only_one_trait_flag_error(self, tmp_path):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            trait_x="body_mass",
            trait_y=None,
            color_by=None,
            plot_output=str(tmp_path / "test.png"),
            json=False,
        )
        svc = Phylomorphospace(args)
        with pytest.raises(PhykitUserError):
            svc.run()

    def test_plot_creates_file(self, default_args, tmp_path):
        default_args.plot_output = str(tmp_path / "test.png")
        svc = Phylomorphospace(default_args)
        svc.run()
        assert Path(default_args.plot_output).exists()
        assert Path(default_args.plot_output).stat().st_size > 0


class TestPlot:
    def test_with_color_by_column(self, tmp_path):
        plot_path = str(tmp_path / "color_col.png")
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            trait_x="body_mass",
            trait_y="brain_size",
            color_by="longevity",
            plot_output=plot_path,
            json=False,
        )
        svc = Phylomorphospace(args)
        svc.run()
        assert Path(plot_path).exists()
        assert Path(plot_path).stat().st_size > 0

    def test_with_color_by_discrete_file(self, tmp_path):
        color_file = tmp_path / "groups.tsv"
        color_file.write_text(
            "raccoon\tgroup_A\n"
            "bear\tgroup_A\n"
            "sea_lion\tgroup_B\n"
            "seal\tgroup_B\n"
            "monkey\tgroup_A\n"
            "cat\tgroup_B\n"
            "weasel\tgroup_A\n"
            "dog\tgroup_B\n"
        )
        plot_path = str(tmp_path / "color_disc.png")
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            trait_x="body_mass",
            trait_y="brain_size",
            color_by=str(color_file),
            plot_output=plot_path,
            json=False,
        )
        svc = Phylomorphospace(args)
        svc.run()
        assert Path(plot_path).exists()
        assert Path(plot_path).stat().st_size > 0


class TestReconstructAncestralScores:
    @pytest.fixture(autouse=True)
    def _make_trait_data(self, tmp_path):
        """Create local trait data so tests don't depend on sample_files."""
        self._trait_file = str(tmp_path / "traits.tsv")
        Path(self._trait_file).write_text(
            "taxon\tbody_mass\tbrain_size\n"
            "raccoon\t1.04\t1.60\n"
            "bear\t2.39\t2.66\n"
            "sea_lion\t2.30\t2.74\n"
            "seal\t1.88\t2.45\n"
            "monkey\t0.60\t1.85\n"
            "cat\t0.56\t1.30\n"
            "weasel\t-0.30\t0.85\n"
            "dog\t1.18\t1.87\n"
        )

    def _make_svc(self):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=self._trait_file,
            trait_x="body_mass",
            trait_y="brain_size",
            color_by=None,
            plot_output="unused.png",
            json=False,
        )
        return Phylomorphospace(args)

    def _get_data(self, svc):
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(
            svc.trait_data_path, tree_tips
        )
        ordered_names = sorted(traits.keys())
        x_idx = trait_names.index("body_mass")
        y_idx = trait_names.index("brain_size")
        data = np.array(
            [[traits[name][x_idx], traits[name][y_idx]] for name in ordered_names]
        )
        return tree, ordered_names, data

    def test_all_nodes_scored(self):
        svc = self._make_svc()
        tree, ordered_names, data = self._get_data(svc)
        node_estimates, node_distances, tree_pruned = (
            svc._reconstruct_ancestral_scores(tree, data, ordered_names)
        )

        # All clades should have estimates
        all_clades = list(tree_pruned.find_clades())
        for clade in all_clades:
            assert id(clade) in node_estimates

    def test_tips_match_input(self):
        svc = self._make_svc()
        tree, ordered_names, data = self._get_data(svc)
        node_estimates, _, tree_pruned = svc._reconstruct_ancestral_scores(
            tree, data, ordered_names
        )

        # Check that tip estimates match input data
        for tip in tree_pruned.get_terminals():
            if tip.name in ordered_names:
                idx = ordered_names.index(tip.name)
                est = node_estimates[id(tip)]
                np.testing.assert_array_almost_equal(est, data[idx])

    def test_root_in_range(self):
        svc = self._make_svc()
        tree, ordered_names, data = self._get_data(svc)
        node_estimates, _, tree_pruned = svc._reconstruct_ancestral_scores(
            tree, data, ordered_names
        )

        root_est = node_estimates[id(tree_pruned.root)]
        # Root estimate should be within the range of tip values
        for dim in range(2):
            assert root_est[dim] >= data[:, dim].min() - 1.0
            assert root_est[dim] <= data[:, dim].max() + 1.0
