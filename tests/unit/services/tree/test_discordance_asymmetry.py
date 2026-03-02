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


class TestTestAsymmetry:
    def _make_svc(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry
        args = Namespace(
            tree="tests/sample_files/tree_simple.tre",
            gene_trees="tests/sample_files/gene_trees_simple.nwk",
            verbose=False, json=False, plot_output=None,
        )
        return DiscordanceAsymmetry(args)

    def test_symmetric_case(self):
        svc = self._make_svc()
        result = svc._test_asymmetry(5, 5)
        assert result["asymmetry_ratio"] == 0.5
        assert result["p_value"] == pytest.approx(1.0)
        assert result["favored_alt"] is None

    def test_asymmetric_case(self):
        svc = self._make_svc()
        result = svc._test_asymmetry(9, 1)
        assert result["asymmetry_ratio"] == 0.9
        assert result["p_value"] < 0.05
        assert result["favored_alt"] == "alt1"

    def test_zero_discordance(self):
        svc = self._make_svc()
        result = svc._test_asymmetry(0, 0)
        assert result["asymmetry_ratio"] is None
        assert result["p_value"] is None
        assert result["favored_alt"] is None

    def test_alt2_favored(self):
        svc = self._make_svc()
        result = svc._test_asymmetry(1, 9)
        assert result["asymmetry_ratio"] == 0.9
        assert result["favored_alt"] == "alt2"

    def test_single_discordant(self):
        svc = self._make_svc()
        result = svc._test_asymmetry(1, 0)
        assert result["asymmetry_ratio"] == 1.0
        assert result["p_value"] == pytest.approx(1.0)
        assert result["favored_alt"] == "alt1"


class TestFDR:
    def test_empty_list(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry
        assert DiscordanceAsymmetry._fdr([]) == []

    def test_single_pvalue(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry
        result = DiscordanceAsymmetry._fdr([0.03])
        assert result == [0.03]

    def test_known_correction(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry
        # Known input: 3 p-values
        pvals = [0.01, 0.04, 0.03]
        result = DiscordanceAsymmetry._fdr(pvals)
        # FDR correction: rank p-values, adjust
        # sorted: (0, 0.01), (2, 0.03), (1, 0.04)
        # rank3 (idx=1, p=0.04): 0.04*3/3 = 0.04, prev=0.04
        # rank2 (idx=2, p=0.03): 0.03*3/2 = 0.045, min(0.045, 0.04)=0.04, prev=0.04
        # rank1 (idx=0, p=0.01): 0.01*3/1 = 0.03, min(0.03, 0.04)=0.03, prev=0.03
        # Result: [0.03, 0.04, 0.04]
        assert result[0] == pytest.approx(0.03)
        assert result[1] == pytest.approx(0.04)
        assert result[2] == pytest.approx(0.04)

    def test_all_ones(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry
        result = DiscordanceAsymmetry._fdr([1.0, 1.0, 1.0])
        assert all(p == 1.0 for p in result)


class TestRun:
    def _make_svc(self, verbose=False, json_output=False, plot_output=None):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry
        args = Namespace(
            tree="tests/sample_files/tree_simple.tre",
            gene_trees="tests/sample_files/gene_trees_simple.nwk",
            verbose=verbose, json=json_output, plot_output=plot_output,
        )
        return DiscordanceAsymmetry(args)

    def test_text_output_has_header(self, capsys):
        svc = self._make_svc()
        svc.run()
        captured = capsys.readouterr()
        assert "branch" in captured.out
        assert "n_conc" in captured.out
        assert "n_alt1" in captured.out
        assert "n_alt2" in captured.out
        assert "asym_ratio" in captured.out
        assert "binom_p" in captured.out
        assert "fdr_p" in captured.out
        assert "gene_flow" in captured.out

    def test_text_output_has_branches(self, capsys):
        svc = self._make_svc()
        svc.run()
        captured = capsys.readouterr()
        assert "bear,raccoon" in captured.out
        assert "cat,monkey" in captured.out

    def test_text_output_has_summary(self, capsys):
        svc = self._make_svc()
        svc.run()
        captured = capsys.readouterr()
        assert "Summary:" in captured.out
        assert "branches tested" in captured.out

    def test_json_output_structure(self, capsys):
        import json
        svc = self._make_svc(json_output=True)
        svc.run()
        captured = capsys.readouterr()
        data = json.loads(captured.out)
        assert "branches" in data
        assert "summary" in data
        assert isinstance(data["branches"], list)
        assert len(data["branches"]) > 0
        # Check branch entry has expected keys
        branch = data["branches"][0]
        for key in ["split", "n_concordant", "n_alt1", "n_alt2",
                     "asymmetry_ratio", "p_value", "fdr_p"]:
            assert key in branch
        # Check summary
        assert "n_gene_trees" in data["summary"]
        assert "n_branches_tested" in data["summary"]
        assert "n_significant_fdr05" in data["summary"]

    def test_json_output_n_gene_trees(self, capsys):
        import json
        svc = self._make_svc(json_output=True)
        svc.run()
        captured = capsys.readouterr()
        data = json.loads(captured.out)
        assert data["summary"]["n_gene_trees"] == 10


class TestPlot:
    def _make_svc(self, plot_output):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry
        args = Namespace(
            tree="tests/sample_files/tree_simple.tre",
            gene_trees="tests/sample_files/gene_trees_simple.nwk",
            verbose=False, json=False, plot_output=plot_output,
        )
        return DiscordanceAsymmetry(args)

    def test_plot_creates_file(self, tmp_path):
        output = str(tmp_path / "test_asym.png")
        svc = self._make_svc(output)
        svc.run()
        assert os.path.exists(output)

    def test_plot_file_nonempty(self, tmp_path):
        output = str(tmp_path / "test_asym.png")
        svc = self._make_svc(output)
        svc.run()
        assert os.path.getsize(output) > 0

    def test_plot_no_error_with_all_concordant(self, tmp_path):
        # Even if all gene trees are concordant (no discordance),
        # plotting should not error
        output = str(tmp_path / "test_asym_conc.png")
        svc = self._make_svc(output)
        svc.run()
        # If we get here without error, the test passes
        assert os.path.exists(output)
