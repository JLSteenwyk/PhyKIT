"""
Unit tests for independent_contrasts (Felsenstein's PIC).

Computes phylogenetically independent contrasts for continuous traits.
Cross-validated against R's ape::pic(). Individual contrasts at
polytomy-resolved nodes may differ due to resolution order, but the
sum of squared contrasts is invariant (verified to match R within
floating-point tolerance).

See tests/r_validation/validate_pic.R for the R validation script.
"""
import pytest
import numpy as np
from argparse import Namespace
from pathlib import Path

from phykit.services.tree.independent_contrasts import IndependentContrasts


@pytest.fixture
def args():
    return Namespace(
        tree="tests/sample_files/tree_simple.tre",
        trait_data="tests/sample_files/tree_simple_traits.tsv",
    )


class TestIndependentContrastsInit:
    def test_init_sets_fields(self, args):
        ic = IndependentContrasts(args)
        assert ic.tree_file_path == args.tree
        assert ic.trait_data_path == args.trait_data
        assert ic.json_output is False

    def test_process_args_json(self):
        args = Namespace(
            tree="t.tre", trait_data="d.tsv", json=True
        )
        ic = IndependentContrasts(args)
        assert ic.json_output is True


class TestPICComputation:
    def test_correct_number_of_contrasts(self, args):
        """n tips should produce n-1 contrasts."""
        ic = IndependentContrasts(args)
        tree = ic.read_tree_file()
        import copy
        tree = copy.deepcopy(tree)
        tree_tips = [t.name for t in tree.get_terminals()]
        tip_traits = ic._parse_trait_data(args.trait_data, tree_tips)
        shared = set(tip_traits.keys())
        tips_to_prune = [t for t in tree_tips if t not in shared]
        if tips_to_prune:
            tree = ic.prune_tree_using_taxa_list(tree, tips_to_prune)
        ic._resolve_polytomies(tree)
        contrasts, _ = ic._compute_pic(tree, tip_traits)
        assert len(contrasts) == len(tip_traits) - 1

    def test_sum_of_squared_contrasts_matches_r(self, args):
        """Sum of squared contrasts should match R's ape::pic() exactly.

        R value: sum(pic(trait_vec, multi2di(tree))^2) = 0.307253
        This is invariant to polytomy resolution order.
        """
        ic = IndependentContrasts(args)
        tree = ic.read_tree_file()
        import copy
        tree = copy.deepcopy(tree)
        tree_tips = [t.name for t in tree.get_terminals()]
        tip_traits = ic._parse_trait_data(args.trait_data, tree_tips)
        shared = set(tip_traits.keys())
        tips_to_prune = [t for t in tree_tips if t not in shared]
        if tips_to_prune:
            tree = ic.prune_tree_using_taxa_list(tree, tips_to_prune)
        ic._resolve_polytomies(tree)
        contrasts, _ = ic._compute_pic(tree, tip_traits)

        ss = sum(c ** 2 for c in contrasts)
        # Cross-validated against R: sum(pic(...)^2) = 0.307253
        assert ss == pytest.approx(0.307253, abs=0.001)

    def test_known_contrasts_from_resolved_clades(self, args):
        """Contrasts from fully resolved clades should match R exactly."""
        ic = IndependentContrasts(args)
        tree = ic.read_tree_file()
        import copy
        tree = copy.deepcopy(tree)
        tree_tips = [t.name for t in tree.get_terminals()]
        tip_traits = ic._parse_trait_data(args.trait_data, tree_tips)
        shared = set(tip_traits.keys())
        tips_to_prune = [t for t in tree_tips if t not in shared]
        if tips_to_prune:
            tree = ic.prune_tree_using_taxa_list(tree, tips_to_prune)
        ic._resolve_polytomies(tree)
        contrasts, node_labels = ic._compute_pic(tree, tip_traits)

        # Find contrasts by their tip composition
        contrast_by_tips = {}
        for c, tips in zip(contrasts, node_labels):
            contrast_by_tips[frozenset(tips)] = c

        # bear-raccoon: R gives -0.264757
        br = frozenset(["bear", "raccoon"])
        assert contrast_by_tips[br] == pytest.approx(-0.264757, abs=0.0001)

        # sea_lion-seal: R gives 0.085732
        ss = frozenset(["sea_lion", "seal"])
        assert contrast_by_tips[ss] == pytest.approx(0.085732, abs=0.0001)

        # monkey-cat: R gives 0.003288
        mc = frozenset(["monkey", "cat"])
        assert contrast_by_tips[mc] == pytest.approx(0.003288, abs=0.0001)

    def test_identical_traits_produce_zero_contrasts(self):
        """If all taxa have the same trait value, all contrasts should be ~0."""
        from Bio import Phylo
        from io import StringIO
        import copy

        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        args = Namespace(tree="dummy", trait_data="dummy")
        ic = IndependentContrasts(args)
        ic._resolve_polytomies(tree)
        tip_traits = {"A": 5.0, "B": 5.0, "C": 5.0, "D": 5.0}
        contrasts, _ = ic._compute_pic(tree, tip_traits)
        for c in contrasts:
            assert abs(c) < 1e-10


class TestPICRun:
    def test_run_text_output(self, args, capsys):
        ic = IndependentContrasts(args)
        ic.run()
        captured = capsys.readouterr()
        assert "Number of contrasts: 7" in captured.out
        assert "Mean absolute contrast:" in captured.out
        assert "Variance of contrasts:" in captured.out

    def test_run_json_output(self, mocker, args):
        args.json = True
        ic = IndependentContrasts(args)
        mocked_json = mocker.patch(
            "phykit.services.tree.independent_contrasts.print_json"
        )
        ic.run()
        payload = mocked_json.call_args.args[0]
        assert payload["n_taxa"] == 8
        assert payload["n_contrasts"] == 7
        assert len(payload["contrasts"]) == 7
        # Sum of squared contrasts matches R
        ss = sum(c["contrast"] ** 2 for c in payload["contrasts"])
        assert ss == pytest.approx(0.307253, abs=0.001)
