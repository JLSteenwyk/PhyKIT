import os
import sys
import tempfile

import pytest
from mock import patch
from argparse import Namespace

from phykit.errors import PhykitUserError

TREE_SIMPLE = "tests/sample_files/tree_simple.tre"
GENE_TREES = "tests/sample_files/gene_trees_simple.nwk"


class TestProcessArgs:
    def test_default_args(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            verbose=False,
            json=False,
            plot_output=None,
        )
        svc = DiscordanceAsymmetry(args)
        assert svc.gene_trees_path == GENE_TREES
        assert svc.verbose is False
        assert svc.json_output is False
        assert svc.plot_output is None

    def test_verbose_arg(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            verbose=True,
            json=False,
            plot_output=None,
        )
        svc = DiscordanceAsymmetry(args)
        assert svc.verbose is True

    def test_json_arg(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            verbose=False,
            json=True,
            plot_output=None,
        )
        svc = DiscordanceAsymmetry(args)
        assert svc.json_output is True

    def test_plot_output_arg(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            verbose=False,
            json=False,
            plot_output="/tmp/test.png",
        )
        svc = DiscordanceAsymmetry(args)
        assert svc.plot_output == "/tmp/test.png"


class TestParseGeneTrees:
    def test_parses_multi_newick(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=False, json=False, plot_output=None,
        )
        svc = DiscordanceAsymmetry(args)
        trees = svc._parse_gene_trees(GENE_TREES)
        assert len(trees) == 10
        for gt in trees:
            tips = [t.name for t in gt.get_terminals()]
            assert len(tips) == 8

    def test_missing_file_raises_error(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        args = Namespace(
            tree=TREE_SIMPLE, gene_trees="nonexistent.nwk",
            verbose=False, json=False, plot_output=None,
        )
        svc = DiscordanceAsymmetry(args)
        with pytest.raises(PhykitUserError):
            svc._parse_gene_trees("nonexistent.nwk")


class TestCountTopologies:
    @pytest.fixture
    def svc(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=False, json=False, plot_output=None,
        )
        return DiscordanceAsymmetry(args)

    def test_counts_with_sample_data(self, svc):
        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        result = svc._count_topologies(species_tree, gene_trees)
        assert isinstance(result, dict)
        assert len(result) > 0
        for branch_key, data in result.items():
            assert "split" in data
            assert "n_concordant" in data
            assert "n_alt1" in data
            assert "n_alt2" in data
            total = data["n_concordant"] + data["n_alt1"] + data["n_alt2"]
            assert total <= 10

    def test_all_branches_present(self, svc):
        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        result = svc._count_topologies(species_tree, gene_trees)
        # Expected branches for the 8-taxon tree:
        # ((raccoon,bear),((sea_lion,seal),((monkey,cat),weasel)),dog)
        # Internal branches produce these canonical splits:
        expected_branches = {
            "bear,raccoon",
            "cat,monkey",
            "cat,monkey,weasel",
            "sea_lion,seal",
            "bear,dog,raccoon",
        }
        assert set(result.keys()) == expected_branches

    def test_concordance_proportions_cross_validate(self, svc):
        """For each branch, gCF = n_concordant / (n_concordant + n_alt1 + n_alt2)
        should be > 0.5 for all branches (majority are concordant in sample data)."""
        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        result = svc._count_topologies(species_tree, gene_trees)
        for branch_key, data in result.items():
            total = data["n_concordant"] + data["n_alt1"] + data["n_alt2"]
            if total == 0:
                continue
            gcf = data["n_concordant"] / total
            assert gcf > 0.5, (
                f"Expected gCF > 0.5 for branch {branch_key}, "
                f"got {gcf:.3f} (conc={data['n_concordant']}, "
                f"alt1={data['n_alt1']}, alt2={data['n_alt2']})"
            )
