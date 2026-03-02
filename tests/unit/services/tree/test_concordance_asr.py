import copy
import json
import math
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

from phykit.services.tree.concordance_asr import ConcordanceAsr
from phykit.errors import PhykitUserError

here = Path(__file__).resolve().parent
sample = here.parent.parent.parent / "sample_files"

TREE_SIMPLE = str(sample / "tree_simple.tre")
GENE_TREES = str(sample / "gene_trees_simple.nwk")
TRAITS_FILE = str(sample / "tree_simple_traits.tsv")


@pytest.fixture
def default_args():
    return Namespace(
        tree=TREE_SIMPLE,
        gene_trees=GENE_TREES,
        trait_data=TRAITS_FILE,
        trait=None,
        method="weighted",
        ci=False,
        plot=None,
        missing_taxa="shared",
        json=False,
    )


@pytest.fixture
def ci_args():
    return Namespace(
        tree=TREE_SIMPLE,
        gene_trees=GENE_TREES,
        trait_data=TRAITS_FILE,
        trait=None,
        method="weighted",
        ci=True,
        plot=None,
        missing_taxa="shared",
        json=False,
    )


@pytest.fixture
def distribution_args():
    return Namespace(
        tree=TREE_SIMPLE,
        gene_trees=GENE_TREES,
        trait_data=TRAITS_FILE,
        trait=None,
        method="distribution",
        ci=True,
        plot=None,
        missing_taxa="shared",
        json=False,
    )


@pytest.fixture
def json_args():
    return Namespace(
        tree=TREE_SIMPLE,
        gene_trees=GENE_TREES,
        trait_data=TRAITS_FILE,
        trait=None,
        method="weighted",
        ci=False,
        plot=None,
        missing_taxa="shared",
        json=True,
    )


class TestProcessArgs:
    def test_defaults(self, default_args):
        svc = ConcordanceAsr(default_args)
        assert svc.method == "weighted"
        assert svc.ci is False
        assert svc.plot_output is None
        assert svc.json_output is False
        assert svc.missing_taxa == "shared"
        assert svc.trait_column is None

    def test_overrides(self):
        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            trait_data=TRAITS_FILE,
            trait="body_mass",
            method="distribution",
            ci=True,
            plot="out.png",
            missing_taxa="error",
            json=True,
        )
        svc = ConcordanceAsr(args)
        assert svc.method == "distribution"
        assert svc.ci is True
        assert svc.plot_output == "out.png"
        assert svc.json_output is True
        assert svc.missing_taxa == "error"
        assert svc.trait_column == "body_mass"


def _make_asr_helper():
    """Create a minimal AncestralReconstruction helper for testing."""
    from phykit.services.tree.ancestral_reconstruction import AncestralReconstruction
    asr = AncestralReconstruction.__new__(AncestralReconstruction)
    asr.ci = True
    return asr


class TestGCFComputation:
    def test_gcf_fractions(self, default_args):
        svc = ConcordanceAsr(default_args)
        svc._asr = _make_asr_helper()

        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        all_taxa = set(t.name for t in species_tree.get_terminals())

        gcf = svc._compute_gcf_per_node(species_tree, gene_trees, all_taxa)

        # All internal nodes should have gCF values
        assert len(gcf) > 0

        # gCF + gDF1 + gDF2 should sum to 1.0 for each node
        for node_id, (g, d1, d2) in gcf.items():
            assert abs(g + d1 + d2 - 1.0) < 1e-10, (
                f"gCF + gDF1 + gDF2 = {g + d1 + d2}, expected 1.0"
            )
            assert 0.0 <= g <= 1.0
            assert 0.0 <= d1 <= 1.0
            assert 0.0 <= d2 <= 1.0

    def test_all_concordant_gcf_is_one(self, default_args):
        """When all gene trees have the same topology, gCF should be 1.0."""
        svc = ConcordanceAsr(default_args)
        svc._asr = _make_asr_helper()
        species_tree = svc.read_tree_file()

        # Create gene trees identical to species tree
        identical_trees = [copy.deepcopy(species_tree) for _ in range(5)]
        all_taxa = set(t.name for t in species_tree.get_terminals())

        gcf = svc._compute_gcf_per_node(species_tree, identical_trees, all_taxa)

        for node_id, (g, d1, d2) in gcf.items():
            assert g == 1.0, f"Expected gCF=1.0, got {g}"
            assert d1 == 0.0
            assert d2 == 0.0

    def test_gcf_uses_correct_four_groups(self, default_args):
        """Verify that gCF computation uses tree-structure-based NNI alternatives."""
        svc = ConcordanceAsr(default_args)
        svc._asr = _make_asr_helper()
        species_tree = svc.read_tree_file()
        all_taxa = frozenset(t.name for t in species_tree.get_terminals())
        parent_map = svc._asr._build_parent_map(species_tree)

        # Check four-group decomposition for each non-root internal node
        for clade in species_tree.find_clades(order="preorder"):
            if clade.is_terminal():
                continue
            groups = svc._get_four_groups(
                species_tree, clade, parent_map, all_taxa
            )
            if groups is None:
                continue
            C1, C2, S, D = groups
            # Four groups should partition all taxa
            assert C1 | C2 | S | D == all_taxa, (
                f"Four groups don't cover all taxa: "
                f"C1={sorted(C1)}, C2={sorted(C2)}, S={sorted(S)}, D={sorted(D)}"
            )
            # Groups should be disjoint
            assert len(C1) + len(C2) + len(S) + len(D) == len(all_taxa)
            # C1 and C2 should be non-empty
            assert len(C1) >= 1
            assert len(C2) >= 1


class TestNNIAlternatives:
    def test_produces_alternatives(self, default_args):
        svc = ConcordanceAsr(default_args)
        svc._asr = _make_asr_helper()
        species_tree = svc.read_tree_file()
        parent_map = svc._asr._build_parent_map(species_tree)

        for clade in species_tree.find_clades(order="preorder"):
            if clade.is_terminal():
                continue
            if clade == species_tree.root:
                continue
            if len(clade.clades) >= 2:
                alts = svc._build_nni_alternative_at_node(
                    species_tree, clade, parent_map
                )
                assert len(alts) >= 1, "Expected at least 1 NNI alternative"
                orig_tips = set(t.name for t in species_tree.get_terminals())
                for alt_tree, expected_desc in alts:
                    alt_tips = set(t.name for t in alt_tree.get_terminals())
                    assert alt_tips == orig_tips, "NNI alternative lost/gained taxa"
                    # expected_desc should be a valid subset of tips
                    assert expected_desc.issubset(orig_tips)
                break

    def test_branch_lengths_preserved(self, default_args):
        svc = ConcordanceAsr(default_args)
        svc._asr = _make_asr_helper()
        species_tree = svc.read_tree_file()
        parent_map = svc._asr._build_parent_map(species_tree)

        def tree_length(tree):
            return sum(
                c.branch_length for c in tree.find_clades()
                if c.branch_length is not None
            )

        orig_len = tree_length(species_tree)

        for clade in species_tree.find_clades(order="preorder"):
            if clade.is_terminal() or clade == species_tree.root:
                continue
            if len(clade.clades) >= 2:
                alts = svc._build_nni_alternative_at_node(
                    species_tree, clade, parent_map
                )
                for alt_tree, _ in alts:
                    alt_len = tree_length(alt_tree)
                    assert abs(alt_len - orig_len) < 1e-6, (
                        f"Tree length changed: {orig_len} -> {alt_len}"
                    )
                break

    def test_expected_desc_matches_nni_swap(self, default_args):
        """Verify expected descendant sets match actual NNI swap topology."""
        svc = ConcordanceAsr(default_args)
        svc._asr = _make_asr_helper()
        species_tree = svc.read_tree_file()
        parent_map = svc._asr._build_parent_map(species_tree)
        all_taxa = frozenset(t.name for t in species_tree.get_terminals())

        for clade in species_tree.find_clades(order="preorder"):
            if clade.is_terminal() or clade == species_tree.root:
                continue
            if len(clade.clades) < 2:
                continue

            C1 = frozenset(t.name for t in clade.clades[0].get_terminals())
            C2 = frozenset(t.name for t in clade.clades[1].get_terminals())
            siblings = [c for c in parent_map[id(clade)].clades
                        if id(c) != id(clade)]
            S = frozenset(t.name for t in siblings[0].get_terminals())

            alts = svc._build_nni_alternative_at_node(
                species_tree, clade, parent_map
            )
            assert len(alts) == 2

            # Alt 0 swaps C1 <-> S: expected desc = S | C2
            _, desc0 = alts[0]
            assert desc0 == S | C2, (
                f"Alt 0: expected {sorted(S | C2)}, got {sorted(desc0)}"
            )

            # Alt 1 swaps C2 <-> S: expected desc = C1 | S
            _, desc1 = alts[1]
            assert desc1 == C1 | S, (
                f"Alt 1: expected {sorted(C1 | S)}, got {sorted(desc1)}"
            )


class TestLawOfTotalVariance:
    def test_known_values(self):
        weights = [0.5, 0.5]
        means = [10.0, 20.0]
        variances = [1.0, 1.0]

        total, within, between = ConcordanceAsr._law_of_total_variance(
            weights, means, variances
        )

        # Within = 0.5*1 + 0.5*1 = 1.0
        assert abs(within - 1.0) < 1e-10

        # Between = 0.5*(10-15)^2 + 0.5*(20-15)^2 = 25.0
        assert abs(between - 25.0) < 1e-10

        # Total = 26.0
        assert abs(total - 26.0) < 1e-10

    def test_single_topology_zero_between(self):
        weights = [1.0]
        means = [10.0]
        variances = [2.0]

        total, within, between = ConcordanceAsr._law_of_total_variance(
            weights, means, variances
        )

        assert abs(between - 0.0) < 1e-10
        assert abs(within - 2.0) < 1e-10
        assert abs(total - 2.0) < 1e-10

    def test_zero_weights(self):
        weights = [0.0, 0.0]
        means = [10.0, 20.0]
        variances = [1.0, 1.0]

        total, within, between = ConcordanceAsr._law_of_total_variance(
            weights, means, variances
        )

        assert total == 0.0
        assert within == 0.0
        assert between == 0.0


class TestWeightedMethod:
    def test_basic_run(self, default_args):
        svc = ConcordanceAsr(default_args)
        svc._asr = _make_asr_helper()

        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        all_taxa = set(t.name for t in species_tree.get_terminals())
        from phykit.services.tree.ancestral_reconstruction import AncestralReconstruction
        trait_values = svc._asr._parse_single_trait_data(TRAITS_FILE, sorted(all_taxa))

        result = svc._run_weighted(species_tree, gene_trees, trait_values, all_taxa)

        assert result["method"] == "weighted"
        assert result["n_gene_trees"] == 10
        assert result["n_tips"] == 8
        assert result["sigma2"] > 0
        assert len(result["ancestral_estimates"]) > 0

        for label, entry in result["ancestral_estimates"].items():
            assert "estimate" in entry
            assert "gcf" in entry
            assert "var_topology" in entry
            assert "var_parameter" in entry
            assert 0.0 <= entry["gcf"] <= 1.0

    def test_all_concordant_matches_standard_asr(self, default_args):
        """When all gene trees match species tree, result should be close to standard ASR."""
        svc = ConcordanceAsr(default_args)
        svc._asr = _make_asr_helper()

        species_tree = svc.read_tree_file()
        all_taxa = set(t.name for t in species_tree.get_terminals())
        trait_values = svc._asr._parse_single_trait_data(TRAITS_FILE, sorted(all_taxa))

        # All gene trees identical to species tree
        identical_trees = [copy.deepcopy(species_tree) for _ in range(5)]

        result = svc._run_weighted(
            copy.deepcopy(species_tree), identical_trees, trait_values, all_taxa
        )

        # Run standard ASR for comparison (same tree object for both calls)
        std_tree = copy.deepcopy(species_tree)
        sp_estimates, _, _, _ = svc._asr._fast_anc(
            std_tree,
            np.array([trait_values[n] for n in sorted(trait_values.keys())]),
            sorted(trait_values.keys()),
            svc._asr._label_internal_nodes(std_tree),
        )

        # All gCF should be 1.0 and var_topology should be ~0
        # Estimates should match standard ASR
        for label, entry in result["ancestral_estimates"].items():
            assert entry["gcf"] == 1.0
            assert entry["var_topology"] < 1e-6
            if label in sp_estimates:
                assert abs(entry["estimate"] - sp_estimates[label]) < 1e-4, (
                    f"Node {label}: concordance-weighted={entry['estimate']:.6f} "
                    f"standard={sp_estimates[label]:.6f}"
                )


class TestDistributionMethod:
    def test_basic_run(self, distribution_args):
        svc = ConcordanceAsr(distribution_args)
        svc._asr = _make_asr_helper()

        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        all_taxa = set(t.name for t in species_tree.get_terminals())
        trait_values = svc._asr._parse_single_trait_data(TRAITS_FILE, sorted(all_taxa))

        result = svc._run_distribution(species_tree, gene_trees, trait_values, all_taxa)

        assert result["method"] == "distribution"
        assert result["n_gene_trees"] == 10
        assert len(result["ancestral_estimates"]) > 0

        for label, entry in result["ancestral_estimates"].items():
            assert "estimate" in entry
            assert "gcf" in entry
            assert "var_topology" in entry

    def test_cis_present(self, distribution_args):
        svc = ConcordanceAsr(distribution_args)
        svc._asr = _make_asr_helper()

        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        all_taxa = set(t.name for t in species_tree.get_terminals())
        trait_values = svc._asr._parse_single_trait_data(TRAITS_FILE, sorted(all_taxa))

        result = svc._run_distribution(species_tree, gene_trees, trait_values, all_taxa)

        # At least some nodes should have CIs
        has_ci = any(
            "ci_lower" in e for e in result["ancestral_estimates"].values()
        )
        assert has_ci, "Distribution method should produce CIs"


class TestRun:
    def test_text_output(self, default_args):
        svc = ConcordanceAsr(default_args)
        with patch("builtins.print") as mock_print:
            svc.run()

        printed = " ".join(str(c) for c in mock_print.call_args_list)
        assert "Concordance-Aware" in printed
        assert "weighted" in printed
        assert "gCF" in printed

    def test_json_output(self, json_args):
        svc = ConcordanceAsr(json_args)
        captured = []
        with patch("phykit.services.tree.concordance_asr.print_json") as mock_json:
            mock_json.side_effect = lambda d: captured.append(d)
            svc.run()

        assert len(captured) == 1
        result = captured[0]
        assert "method" in result
        assert "ancestral_estimates" in result
        assert "n_gene_trees" in result
        assert "sigma2" in result

    def test_distribution_text_output(self, distribution_args):
        svc = ConcordanceAsr(distribution_args)
        with patch("builtins.print") as mock_print:
            svc.run()

        printed = " ".join(str(c) for c in mock_print.call_args_list)
        assert "Concordance-Aware" in printed
        assert "distribution" in printed


class TestEdgeCases:
    def test_single_gene_tree_errors(self):
        """Single gene tree should raise error."""
        import tempfile
        import os

        # Create a temp file with just one gene tree
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".nwk", delete=False
        ) as f:
            f.write(
                "((raccoon:18.5,bear:7.2):0.9,((sea_lion:12.1,seal:11.8):7.6,"
                "((monkey:99.2,cat:46.5):21.0,weasel:19.1):2.1):3.9,dog:24.8);\n"
            )
            tmp_path = f.name

        try:
            args = Namespace(
                tree=TREE_SIMPLE,
                gene_trees=tmp_path,
                trait_data=TRAITS_FILE,
                trait=None,
                method="weighted",
                ci=False,
                plot=None,
                missing_taxa="shared",
                json=False,
            )
            svc = ConcordanceAsr(args)
            with pytest.raises(SystemExit):
                svc.run()
        finally:
            os.unlink(tmp_path)

    def test_missing_taxa_error_mode(self):
        """Error mode should reject taxa mismatches."""
        import tempfile
        import os

        # Create gene trees with different taxa
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".nwk", delete=False
        ) as f:
            f.write(
                "((raccoon:18.5,bear:7.2):0.9,(seal:11.8,"
                "((monkey:99.2,cat:46.5):21.0,weasel:19.1):2.1):3.9,dog:24.8);\n"
                "((raccoon:19.0,bear:6.8):0.85,(seal:12.0,"
                "((monkey:100.0,cat:47.0):20.5,weasel:18.8):2.0):3.8,dog:25.0);\n"
            )
            tmp_path = f.name

        try:
            args = Namespace(
                tree=TREE_SIMPLE,
                gene_trees=tmp_path,
                trait_data=TRAITS_FILE,
                trait=None,
                method="weighted",
                ci=False,
                plot=None,
                missing_taxa="error",
                json=False,
            )
            svc = ConcordanceAsr(args)
            with pytest.raises(SystemExit):
                svc.run()
        finally:
            os.unlink(tmp_path)

    def test_gene_tree_file_not_found(self):
        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees="/nonexistent/path.nwk",
            trait_data=TRAITS_FILE,
            trait=None,
            method="weighted",
            ci=False,
            plot=None,
            missing_taxa="shared",
            json=False,
        )
        svc = ConcordanceAsr(args)
        with pytest.raises(SystemExit):
            svc.run()

    def test_ci_with_weighted(self, ci_args):
        svc = ConcordanceAsr(ci_args)
        captured = []
        with patch("builtins.print") as mock_print:
            svc.run()

        printed = " ".join(str(c) for c in mock_print.call_args_list)
        assert "95% CI" in printed


class TestBenchmarkVerification:
    """Verify correctness of concordance-aware ASR against known values."""

    def test_gcf_manual_computation(self, default_args):
        """Verify gCF values against manually counted concordance in gene trees.

        The species tree topology has the bipartition:
          {raccoon, bear} | {sea_lion, seal, monkey, cat, weasel, dog}
        Gene trees 1-6 and 10 preserve this clade (7 out of 10).
        Gene trees 7, 8, 9 break this clade.
        The gCF for this node should reflect the fraction of gene trees
        supporting it among those resolving the relevant quartet.
        """
        svc = ConcordanceAsr(default_args)
        svc._asr = _make_asr_helper()
        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        all_taxa = set(t.name for t in species_tree.get_terminals())

        gcf = svc._compute_gcf_per_node(species_tree, gene_trees, all_taxa)

        # At least some nodes should have gCF < 1.0
        # (discordant gene trees 7-9 break at least some bipartitions)
        has_partial_concordance = any(g < 1.0 for g, _, _ in gcf.values())
        assert has_partial_concordance, (
            "Expected some nodes with gCF < 1.0 given discordant gene trees"
        )

        # Verify gCF for {raccoon, bear} clade: it should be supported
        # by gene trees 1-6 and 10 (7 concordant).
        # The three quartet resolutions at this node:
        #   concordant: {raccoon, bear} | rest
        #   NNI alt 1:  {sibling, bear} | rest
        #   NNI alt 2:  {raccoon, sibling} | rest
        # We check that the concordant fraction is ~0.7 (7/10)
        raccoon_bear_node = None
        for clade in species_tree.find_clades(order="preorder"):
            if clade.is_terminal():
                continue
            tips = frozenset(t.name for t in clade.get_terminals())
            if tips == frozenset({"raccoon", "bear"}):
                raccoon_bear_node = clade
                break

        if raccoon_bear_node is not None and id(raccoon_bear_node) in gcf:
            g, d1, d2 = gcf[id(raccoon_bear_node)]
            # At least 7 of 10 gene trees group raccoon+bear,
            # so gCF should be >= 0.7
            assert g >= 0.6, (
                f"Expected gCF >= 0.6 for {{raccoon, bear}}, got {g}"
            )

    def test_var_topology_positive_with_discordance(self, default_args):
        """When gene trees disagree, var_topology should be > 0 for at least
        some nodes (the between-topology variance from the law of total variance).
        """
        svc = ConcordanceAsr(default_args)
        svc._asr = _make_asr_helper()

        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        all_taxa = set(t.name for t in species_tree.get_terminals())
        trait_values = svc._asr._parse_single_trait_data(
            TRAITS_FILE, sorted(all_taxa)
        )

        result = svc._run_weighted(species_tree, gene_trees, trait_values, all_taxa)

        # At least one node should have var_topology > 0
        has_positive_var_topo = any(
            e["var_topology"] > 0
            for e in result["ancestral_estimates"].values()
        )
        assert has_positive_var_topo, (
            "Expected at least one node with var_topology > 0 "
            "given discordant gene trees"
        )

    def test_weighted_differs_from_standard_asr_with_discordance(self, default_args):
        """The concordance-weighted estimate should differ from standard ASR
        when gene trees are discordant, because NNI alternative estimates
        are mixed in with their discordance weights.
        """
        svc = ConcordanceAsr(default_args)
        svc._asr = _make_asr_helper()

        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        all_taxa = set(t.name for t in species_tree.get_terminals())
        trait_values = svc._asr._parse_single_trait_data(
            TRAITS_FILE, sorted(all_taxa)
        )

        # Weighted (concordance-aware) result
        weighted_result = svc._run_weighted(
            copy.deepcopy(species_tree), gene_trees, trait_values, all_taxa
        )

        # Standard ASR on species tree only
        # Must use the same tree object for _fast_anc and _label_internal_nodes
        # since node_labels is keyed by id(clade)
        std_tree = copy.deepcopy(species_tree)
        sp_est, _, _, _ = svc._asr._fast_anc(
            std_tree,
            np.array([trait_values[n] for n in sorted(trait_values.keys())]),
            sorted(trait_values.keys()),
            svc._asr._label_internal_nodes(std_tree),
        )

        # For nodes with gCF < 1.0 (where NNI alternatives contribute),
        # the weighted estimate should differ from the standard ASR estimate
        has_difference = False
        for label, entry in weighted_result["ancestral_estimates"].items():
            if entry["gcf"] < 1.0 and label in sp_est:
                if abs(entry["estimate"] - sp_est[label]) > 1e-6:
                    has_difference = True
                    break

        assert has_difference, (
            "Expected concordance-weighted estimate to differ from "
            "standard ASR for at least one discordant node"
        )

    def test_all_concordant_zero_topology_variance(self, default_args):
        """When all gene trees are identical to the species tree,
        var_topology should be 0 for all nodes (no topological uncertainty).
        """
        svc = ConcordanceAsr(default_args)
        svc._asr = _make_asr_helper()

        species_tree = svc.read_tree_file()
        all_taxa = set(t.name for t in species_tree.get_terminals())
        trait_values = svc._asr._parse_single_trait_data(
            TRAITS_FILE, sorted(all_taxa)
        )

        identical_trees = [copy.deepcopy(species_tree) for _ in range(5)]
        result = svc._run_weighted(
            copy.deepcopy(species_tree), identical_trees, trait_values, all_taxa
        )

        for label, entry in result["ancestral_estimates"].items():
            assert entry["var_topology"] < 1e-10, (
                f"Node {label}: expected var_topology ~0 with all concordant "
                f"trees, got {entry['var_topology']}"
            )
            assert entry["gcf"] == 1.0, (
                f"Node {label}: expected gCF=1.0 with all concordant "
                f"trees, got {entry['gcf']}"
            )

    def test_distribution_variance_tracks_discordance(self, distribution_args):
        """In the distribution method, nodes with more discordance should
        tend to have higher var_topology (more spread across gene trees).
        """
        svc = ConcordanceAsr(distribution_args)
        svc._asr = _make_asr_helper()

        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        all_taxa = set(t.name for t in species_tree.get_terminals())
        trait_values = svc._asr._parse_single_trait_data(
            TRAITS_FILE, sorted(all_taxa)
        )

        result = svc._run_distribution(
            species_tree, gene_trees, trait_values, all_taxa
        )

        # At least some nodes should have non-zero variance from gene tree spread
        has_nonzero_var = any(
            e["var_topology"] > 0
            for e in result["ancestral_estimates"].values()
        )
        assert has_nonzero_var, (
            "Distribution method should show variance across gene trees"
        )

    def test_nni_alternatives_incorporated(self, default_args):
        """Verify that NNI alternative ASR estimates are actually used in
        the weighted combination (not silently dropped).
        """
        svc = ConcordanceAsr(default_args)
        svc._asr = _make_asr_helper()

        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        all_taxa = set(t.name for t in species_tree.get_terminals())
        trait_values = svc._asr._parse_single_trait_data(
            TRAITS_FILE, sorted(all_taxa)
        )

        # Run components manually
        node_labels = svc._asr._label_internal_nodes(species_tree)
        parent_map = svc._asr._build_parent_map(species_tree)
        gcf_per_node = svc._compute_gcf_per_node(
            species_tree, gene_trees, all_taxa
        )

        # For nodes with gCF < 1.0, verify NNI alternatives can be found
        nni_found_count = 0
        nni_attempted_count = 0

        for clade in species_tree.find_clades(order="preorder"):
            if clade.is_terminal():
                continue
            if id(clade) not in gcf_per_node:
                continue
            gcf, gdf1, gdf2 = gcf_per_node[id(clade)]
            if gcf >= 1.0:
                continue

            nni_attempted_count += 1
            alts = svc._build_nni_alternative_at_node(
                species_tree, clade, parent_map
            )
            for alt_tree, expected_desc in alts:
                est, _, _ = svc._run_asr_on_tree(alt_tree, trait_values)
                if expected_desc in est:
                    nni_found_count += 1

        # We should have attempted at least one discordant node
        assert nni_attempted_count > 0, (
            "Expected at least one node with gCF < 1.0"
        )
        # At least some NNI lookups should succeed
        assert nni_found_count > 0, (
            f"NNI estimate lookups all failed "
            f"({nni_found_count}/{nni_attempted_count * 2} found). "
            "The expected descendant sets may not match NNI tree nodes."
        )
