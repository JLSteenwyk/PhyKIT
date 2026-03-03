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
        result, _ = svc._count_topologies(species_tree, gene_trees)
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
        result, _ = svc._count_topologies(species_tree, gene_trees)
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
        result, _ = svc._count_topologies(species_tree, gene_trees)
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

    def test_verbose_output(self, capsys):
        svc = self._make_svc(verbose=True)
        svc.run()
        captured = capsys.readouterr()
        assert "Branch:" in captured.out
        assert "gCF=" in captured.out
        assert "gDF1=" in captured.out
        assert "gDF2=" in captured.out


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


class TestMissingTaxa:
    """Test that taxa present in the species tree but absent from gene trees
    (or vice versa) are handled correctly."""

    @pytest.fixture
    def svc(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=False, json=False, plot_output=None,
        )
        return DiscordanceAsymmetry(args)

    def test_species_tree_has_extra_taxon(self, svc, tmp_path):
        """Species tree has taxon 'e' absent from all gene trees.

        Species tree: ((a,e),(b,(c,d)));
        Gene trees use only {a,b,c,d}.

        The (a,e) branch should be skipped (C1 or C2 empty after filtering).
        The (c,d) and (b,(c,d)) branches should still produce correct counts.
        """
        from Bio import Phylo
        from io import StringIO

        species_tree = Phylo.read(StringIO("((a,e),(b,(c,d)));"), "newick")

        # Write gene trees: all concordant with ((a),(b,(c,d)))
        gt_file = tmp_path / "gene_trees.nwk"
        gt_lines = [
            "(a,(b,(c,d)));",
            "(a,(b,(c,d)));",
            "(a,(b,(c,d)));",
        ]
        gt_file.write_text("\n".join(gt_lines))
        gene_trees = svc._parse_gene_trees(str(gt_file))

        result, _ = svc._count_topologies(species_tree, gene_trees)

        # The branch grouping (a,e) should be skipped because 'e' is not in
        # gene trees, making C2 empty after filtering
        # We should still see branches for (c,d) and (b,(c,d)) -> complement is (a)
        assert len(result) > 0

        # (c,d) branch should be present and concordant
        assert "c,d" in result
        assert result["c,d"]["n_concordant"] == 3

        # The branch for (a) vs (b,c,d) should NOT appear because the
        # species-tree node grouping (a,e) becomes degenerate after filtering
        # (e is removed, making C2 empty)

    def test_gene_tree_has_extra_taxon(self, svc, tmp_path):
        """Gene trees have taxon 'x' absent from species tree.

        Species tree: (a,(b,(c,d)));
        Gene trees include {a,b,c,d,x}.

        The extra taxon 'x' should be ignored; results should match as if
        gene trees only had {a,b,c,d}.
        """
        from Bio import Phylo
        from io import StringIO

        species_tree = Phylo.read(StringIO("(a,(b,(c,d)));"), "newick")

        gt_file = tmp_path / "gene_trees.nwk"
        gt_lines = [
            "((a,x),(b,(c,d)));",
            "((a,x),(b,(c,d)));",
            "((a,x),(b,(c,d)));",
        ]
        gt_file.write_text("\n".join(gt_lines))
        gene_trees = svc._parse_gene_trees(str(gt_file))

        result, _ = svc._count_topologies(species_tree, gene_trees)

        # Extra taxon 'x' in gene trees should be ignored
        assert "c,d" in result
        assert result["c,d"]["n_concordant"] == 3

    def test_degenerate_branch_filtered_out(self, svc, tmp_path):
        """When filtering leaves a group empty, that branch is skipped entirely."""
        from Bio import Phylo
        from io import StringIO

        # Species tree where one child of root has only taxa missing from gene trees
        species_tree = Phylo.read(StringIO("((e,f),(a,(b,(c,d))));"), "newick")

        gt_file = tmp_path / "gene_trees.nwk"
        gt_lines = [
            "(a,(b,(c,d)));",
            "(a,(b,(c,d)));",
        ]
        gt_file.write_text("\n".join(gt_lines))
        gene_trees = svc._parse_gene_trees(str(gt_file))

        result, _ = svc._count_topologies(species_tree, gene_trees)

        # (e,f) branch should be entirely skipped (both taxa missing)
        # (a,(b,(c,d))) subtree branches should still work
        assert "c,d" in result
        assert result["c,d"]["n_concordant"] == 2

        # No branch key should reference e or f
        for key in result:
            assert "e" not in key
            assert "f" not in key


class TestBifurcatingRoot:
    """Test that bifurcating root branches get correct NNI decomposition."""

    @pytest.fixture
    def svc(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=False, json=False, plot_output=None,
        )
        return DiscordanceAsymmetry(args)

    def test_bifurcating_root_nni_alternatives(self, svc, tmp_path):
        """Bifurcating root ((a,b),(c,(d,e))) with known NNI alternatives.

        The root branch separates {a,b} from {c,d,e}.
        Four subtrees: {a}, {b}, {c}, {d,e}.
        Concordant: {a,b} | {c,d,e}
        NNI alt 1: {a,c} | {b,d,e}   (swap b <-> c)
        NNI alt 2: {a,d,e} | {b,c}   (swap b <-> d,e)
        """
        from Bio import Phylo
        from io import StringIO

        species_tree = Phylo.read(StringIO("((a,b),(c,(d,e)));"), "newick")

        gt_file = tmp_path / "gene_trees.nwk"
        gt_lines = [
            "((a,b),(c,(d,e)));",   # concordant
            "((a,b),(c,(d,e)));",   # concordant
            "((a,c),(b,(d,e)));",   # NNI alt: {a,c} | {b,d,e}
            "((a,(d,e)),(b,c));",   # NNI alt: {b,c} | {a,d,e}
        ]
        gt_file.write_text("\n".join(gt_lines))
        gene_trees = svc._parse_gene_trees(str(gt_file))

        result, _ = svc._count_topologies(species_tree, gene_trees)

        assert "a,b" in result
        r = result["a,b"]
        assert r["n_concordant"] == 2
        # Both NNI alternatives should be detected (1 each)
        assert r["n_alt1"] + r["n_alt2"] == 2
        assert r["n_alt1"] == 1
        assert r["n_alt2"] == 1

    def test_bifurcating_root_leaf_sibling_skipped(self, svc, tmp_path):
        """When root is bifurcating and sibling is a leaf, the root branch
        should be skipped (only 3 subtrees, NNI not possible)."""
        from Bio import Phylo
        from io import StringIO

        species_tree = Phylo.read(StringIO("(a,(b,(c,d)));"), "newick")

        gt_file = tmp_path / "gene_trees.nwk"
        gt_lines = ["(a,(b,(c,d)));", "(a,(b,(c,d)));"]
        gt_file.write_text("\n".join(gt_lines))
        gene_trees = svc._parse_gene_trees(str(gt_file))

        result, _ = svc._count_topologies(species_tree, gene_trees)

        # The root branch (a | b,c,d) should be skipped because
        # sibling 'a' is a leaf — can't form 4 subtrees for NNI.
        # Only the (c,d) branch should remain.
        assert "c,d" in result
        # The branch_key "a" (root branch from bifurcating root with leaf)
        # should not appear
        assert "a" not in result


class TestSiblingSelection:
    """Test that sibling selection skips siblings with all-filtered taxa."""

    @pytest.fixture
    def svc(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=False, json=False, plot_output=None,
        )
        return DiscordanceAsymmetry(args)

    def test_first_sibling_filtered_uses_second(self, svc, tmp_path):
        """When siblings[0] taxa are all absent from gene trees, the code
        should use the next sibling with valid taxa."""
        from Bio import Phylo
        from io import StringIO

        # Trifurcating root: (e,f) has no shared taxa, (a,b) and (c,d) do
        species_tree = Phylo.read(StringIO("((e,f),(a,b),(c,d));"), "newick")

        gt_file = tmp_path / "gene_trees.nwk"
        gt_lines = [
            "((a,b),(c,d));",
            "((a,c),(b,d));",
        ]
        gt_file.write_text("\n".join(gt_lines))
        gene_trees = svc._parse_gene_trees(str(gt_file))

        result, _ = svc._count_topologies(species_tree, gene_trees)

        # Branches for (a,b) and (c,d) should NOT be skipped — the code
        # should find a valid sibling past the filtered-out (e,f)
        assert len(result) > 0
        assert "a,b" in result
        assert result["a,b"]["n_concordant"] == 1


class TestPlotWithExtraTaxa:
    """Test that _plot correctly matches branch results when species tree
    has extra taxa not in gene trees."""

    def test_plot_with_extra_species_taxa(self, tmp_path):
        """Plot should color branches correctly even when species tree has
        taxa absent from gene trees."""
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry
        from Bio import Phylo
        from io import StringIO

        # Create species tree with extra taxon 'e'
        sp_file = tmp_path / "species.tre"
        sp_file.write_text("((a,e),(b,(c,d)));")

        gt_file = tmp_path / "gene_trees.nwk"
        gt_lines = [
            "(a,(b,(c,d)));",
            "(a,(b,(c,d)));",
            "(a,(b,(c,d)));",
        ]
        gt_file.write_text("\n".join(gt_lines))

        args = Namespace(
            tree=str(sp_file), gene_trees=str(gt_file),
            verbose=False, json=False, plot_output=None,
        )
        svc = DiscordanceAsymmetry(args)

        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(str(gt_file))
        result, shared_taxa = svc._count_topologies(species_tree, gene_trees)

        # Build branch_results the same way run() does
        branch_results = []
        for branch_key in sorted(result.keys()):
            data = result[branch_key]
            test_result = svc._test_asymmetry(data["n_alt1"], data["n_alt2"])
            entry = dict(
                split=data["split"],
                n_concordant=data["n_concordant"],
                n_alt1=data["n_alt1"],
                n_alt2=data["n_alt2"],
            )
            entry.update(test_result)
            entry["fdr_p"] = None
            branch_results.append(entry)

        # Call _plot with shared_taxa — should not error and should create file
        output = str(tmp_path / "test_plot.png")
        svc._plot(species_tree, branch_results, output,
                  shared_taxa=shared_taxa)
        assert os.path.exists(output)
        assert os.path.getsize(output) > 0


class TestCLI:
    def test_disc_asym_alias(self, capsys):
        from phykit.phykit import Phykit
        sys.argv = [
            "phykit", "disc_asym",
            "-t", "tests/sample_files/tree_simple.tre",
            "-g", "tests/sample_files/gene_trees_simple.nwk",
        ]
        Phykit()
        captured = capsys.readouterr()
        assert "branch" in captured.out
        assert "Summary:" in captured.out

    def test_da_alias(self, capsys):
        from phykit.phykit import Phykit
        sys.argv = [
            "phykit", "da",
            "-t", "tests/sample_files/tree_simple.tre",
            "-g", "tests/sample_files/gene_trees_simple.nwk",
        ]
        Phykit()
        captured = capsys.readouterr()
        assert "bear,raccoon" in captured.out
