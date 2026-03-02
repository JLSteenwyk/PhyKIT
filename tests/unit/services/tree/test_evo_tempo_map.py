import os
import sys
import tempfile

import numpy as np
import pytest
from mock import patch
from argparse import Namespace

from phykit.errors import PhykitUserError

TREE_SIMPLE = "tests/sample_files/tree_simple.tre"
GENE_TREES = "tests/sample_files/gene_trees_simple.nwk"


class TestProcessArgs:
    def test_default_args(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            verbose=False,
            json=False,
            plot_output=None,
        )
        svc = EvoTempoMap(args)
        assert svc.gene_trees_path == GENE_TREES
        assert svc.verbose is False
        assert svc.json_output is False
        assert svc.plot_output is None

    def test_verbose_arg(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            verbose=True,
            json=False,
            plot_output=None,
        )
        svc = EvoTempoMap(args)
        assert svc.verbose is True

    def test_json_arg(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            verbose=False,
            json=True,
            plot_output=None,
        )
        svc = EvoTempoMap(args)
        assert svc.json_output is True

    def test_plot_output_arg(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            verbose=False,
            json=False,
            plot_output="/tmp/test.png",
        )
        svc = EvoTempoMap(args)
        assert svc.plot_output == "/tmp/test.png"


class TestParseGeneTrees:
    def test_parses_multi_newick(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=False, json=False, plot_output=None,
        )
        svc = EvoTempoMap(args)
        trees = svc._parse_gene_trees(GENE_TREES)
        assert len(trees) == 10
        for gt in trees:
            tips = [t.name for t in gt.get_terminals()]
            assert len(tips) == 8

    def test_file_not_found_error(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        args = Namespace(
            tree=TREE_SIMPLE, gene_trees="nonexistent.nwk",
            verbose=False, json=False, plot_output=None,
        )
        svc = EvoTempoMap(args)
        with pytest.raises(PhykitUserError):
            svc._parse_gene_trees("nonexistent.nwk")

    def test_validates_branch_lengths(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        no_bl = "((A,B),(C,D));\n((E,F),(G,H));\n"
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".nwk", delete=False
        ) as f:
            f.write(no_bl)
            tmp_path = f.name
        try:
            args = Namespace(
                tree=TREE_SIMPLE, gene_trees=tmp_path,
                verbose=False, json=False, plot_output=None,
            )
            svc = EvoTempoMap(args)
            with pytest.raises(PhykitUserError) as excinfo:
                svc._parse_and_validate_gene_trees()
            assert "branch lengths" in excinfo.value.messages[1]
        finally:
            os.unlink(tmp_path)

    def test_requires_at_least_2_gene_trees(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        single = "((raccoon:19.2,bear:6.8):0.85,((sea_lion:12.0,seal:12.0):7.5,((monkey:100.9,cat:47.1):20.6,weasel:18.9):2.1):3.9,dog:25.5);\n"
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".nwk", delete=False
        ) as f:
            f.write(single)
            tmp_path = f.name
        try:
            args = Namespace(
                tree=TREE_SIMPLE, gene_trees=tmp_path,
                verbose=False, json=False, plot_output=None,
            )
            svc = EvoTempoMap(args)
            with pytest.raises(PhykitUserError) as excinfo:
                svc._parse_and_validate_gene_trees()
            assert "At least 2" in excinfo.value.messages[0]
        finally:
            os.unlink(tmp_path)


class TestClassifyGeneTrees:
    @pytest.fixture
    def svc(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap
        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=False, json=False, plot_output=None,
        )
        return EvoTempoMap(args)

    def test_classify_returns_dict_per_branch(self, svc):
        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        result = svc._classify_gene_trees(species_tree, gene_trees)
        assert len(result) > 0
        for branch_label, data in result.items():
            assert "concordant_lengths" in data
            assert "discordant_lengths" in data
            assert "n_concordant" in data
            assert "n_discordant" in data
            assert "split" in data

    def test_concordant_plus_discordant_leq_total(self, svc):
        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        result = svc._classify_gene_trees(species_tree, gene_trees)
        n_gene_trees = len(gene_trees)
        for branch_label, data in result.items():
            total = data["n_concordant"] + data["n_discordant"]
            assert total <= n_gene_trees

    def test_all_concordant_when_identical_trees(self, svc):
        species_tree = svc.read_tree_file()
        identical_trees = [svc.read_tree_file() for _ in range(3)]
        result = svc._classify_gene_trees(species_tree, identical_trees)
        for branch_label, data in result.items():
            assert data["n_concordant"] == 3
            assert data["n_discordant"] == 0

    def test_branch_lengths_are_positive(self, svc):
        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        result = svc._classify_gene_trees(species_tree, gene_trees)
        for branch_label, data in result.items():
            for bl in data["concordant_lengths"]:
                assert bl > 0
            for bl in data["discordant_lengths"]:
                assert bl > 0

    def test_cross_check_with_concordance_asr_gcf(self, svc):
        """The gCF (fraction concordant) at each branch should match
        concordance_asr's computation."""
        from phykit.services.tree.concordance_asr import ConcordanceAsr
        from phykit.services.tree.ancestral_reconstruction import AncestralReconstruction

        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        etm_result = svc._classify_gene_trees(species_tree, gene_trees)

        # Get gCF from concordance_asr
        casr_args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            trait_data="tests/sample_files/tree_simple_traits.tsv",
            trait=None, method="weighted", ci=False, plot=None,
            missing_taxa="shared", json=False,
        )
        casr = ConcordanceAsr(casr_args)
        # Initialize the _asr helper that _compute_gcf_per_node depends on
        casr._asr = AncestralReconstruction.__new__(AncestralReconstruction)
        casr._asr.ci = True
        casr_species_tree = casr.read_tree_file()
        casr_gene_trees = casr._parse_gene_trees(GENE_TREES)
        all_taxa = sorted(set(
            t.name for t in casr_species_tree.get_terminals()
        ))
        gcf_result = casr._compute_gcf_per_node(
            casr_species_tree, casr_gene_trees, all_taxa
        )

        # For each species tree branch, compute gCF from ETM and compare
        for branch_key, etm_data in etm_result.items():
            n_conc = etm_data["n_concordant"]
            n_disc = etm_data["n_discordant"]
            total = n_conc + n_disc
            if total == 0:
                continue
            etm_gcf = n_conc / total

            # Find matching node in concordance_asr by descendant set
            split_set = frozenset(etm_data["split"])
            matched = False
            for node_id, (gcf, gdf1, gdf2) in gcf_result.items():
                for clade in casr_species_tree.find_clades():
                    if id(clade) == node_id:
                        node_tips = frozenset(
                            t.name for t in clade.get_terminals()
                        )
                        # canonical split might be the complement
                        complement = frozenset(all_taxa) - node_tips
                        canon = min(node_tips, complement, key=len) if len(node_tips) != len(complement) else min(node_tips, complement, key=lambda s: sorted(s))
                        if canon == split_set:
                            assert abs(etm_gcf - gcf) < 0.01, (
                                f"gCF mismatch at {split_set}: "
                                f"ETM={etm_gcf:.3f}, concordance_asr={gcf:.3f}"
                            )
                            matched = True
                            break
                if matched:
                    break


class TestStatisticalTests:
    @pytest.fixture
    def svc(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap
        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=False, json=False, plot_output=None,
        )
        return EvoTempoMap(args)

    def test_mann_whitney_returns_U_and_pvalue(self, svc):
        conc = [10.0, 12.0, 11.0, 13.0, 10.5]
        disc = [5.0, 6.0, 4.5, 5.5]
        result = svc._test_branch(conc, disc)
        assert "mann_whitney_U" in result
        assert "mann_whitney_p" in result
        assert "permutation_p" in result
        assert result["mann_whitney_p"] < 0.05

    def test_permutation_p_between_0_and_1(self, svc):
        conc = [10.0, 12.0, 11.0, 13.0, 10.5]
        disc = [5.0, 6.0, 4.5, 5.5]
        result = svc._test_branch(conc, disc)
        assert 0.0 <= result["permutation_p"] <= 1.0

    def test_no_difference_gives_high_p(self, svc):
        rng = np.random.default_rng(42)
        vals = rng.normal(10, 1, 20).tolist()
        conc = vals[:10]
        disc = vals[10:]
        result = svc._test_branch(conc, disc)
        assert result["mann_whitney_p"] > 0.05

    def test_insufficient_data_returns_none(self, svc):
        conc = [10.0, 12.0, 11.0]
        disc = [5.0]
        result = svc._test_branch(conc, disc)
        assert result["mann_whitney_p"] is None
        assert result["permutation_p"] is None

    def test_empty_discordant_returns_none(self, svc):
        conc = [10.0, 12.0, 11.0]
        disc = []
        result = svc._test_branch(conc, disc)
        assert result["mann_whitney_p"] is None
        assert result["permutation_p"] is None

    def test_summary_stats_correct(self, svc):
        conc = [10.0, 20.0, 30.0]
        disc = [5.0, 15.0, 25.0]
        result = svc._test_branch(conc, disc)
        assert abs(result["concordant_mean"] - 20.0) < 1e-10
        assert abs(result["concordant_median"] - 20.0) < 1e-10
        assert abs(result["discordant_mean"] - 15.0) < 1e-10
        assert abs(result["discordant_median"] - 15.0) < 1e-10


class TestFDRCorrection:
    def test_fdr_correction(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap
        p_values = [0.01, 0.04, 0.03, 0.20]
        corrected = EvoTempoMap._fdr(p_values)
        assert len(corrected) == 4
        for orig, corr in zip(p_values, corrected):
            assert corr >= orig
        for c in corrected:
            assert c <= 1.0

    def test_fdr_empty_list(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap
        assert EvoTempoMap._fdr([]) == []

    def test_fdr_single_value(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap
        assert EvoTempoMap._fdr([0.05]) == [0.05]

    def test_fdr_matches_relative_rate_test(self):
        """Cross-check: our FDR should match relative_rate_test._fdr."""
        from phykit.services.tree.evo_tempo_map import EvoTempoMap
        from phykit.services.tree.relative_rate_test import RelativeRateTest
        p_values = [0.001, 0.01, 0.05, 0.10, 0.50]
        etm_corrected = EvoTempoMap._fdr(p_values)
        rrt_corrected = RelativeRateTest._fdr(p_values)
        for a, b in zip(etm_corrected, rrt_corrected):
            assert abs(a - b) < 1e-10


class TestGlobalTreeness:
    @pytest.fixture
    def svc(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap
        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=False, json=False, plot_output=None,
        )
        return EvoTempoMap(args)

    def test_treeness_matches_base_class(self, svc):
        """Cross-check: treeness computed here should match
        Tree.calculate_treeness() for each gene tree."""
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        for gt in gene_trees:
            expected = svc.calculate_treeness(gt)
            actual = svc._compute_treeness(gt)
            assert abs(actual - expected) < 1e-10, (
                f"Treeness mismatch: expected {expected}, got {actual}"
            )

    def test_global_treeness_returns_structure(self, svc):
        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        result = svc._compute_global_treeness(species_tree, gene_trees)
        assert "treeness_concordant" in result
        assert "treeness_discordant" in result
        assert "treeness_U_p" in result
        assert "n" in result["treeness_concordant"]
        assert "mean" in result["treeness_concordant"]

    def test_all_concordant_treeness(self, svc):
        species_tree = svc.read_tree_file()
        identical_trees = [svc.read_tree_file() for _ in range(3)]
        result = svc._compute_global_treeness(species_tree, identical_trees)
        assert result["treeness_concordant"]["n"] == 3
        assert result["treeness_discordant"]["n"] == 0
        assert result["treeness_U_p"] is None  # Can't test with 0 discordant

    def test_concordant_plus_discordant_equals_total(self, svc):
        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        result = svc._compute_global_treeness(species_tree, gene_trees)
        total = result["treeness_concordant"]["n"] + result["treeness_discordant"]["n"]
        assert total == 10  # 10 gene trees in the sample file

    def test_treeness_values_between_0_and_1(self, svc):
        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        result = svc._compute_global_treeness(species_tree, gene_trees)
        for group in ["treeness_concordant", "treeness_discordant"]:
            if result[group]["mean"] is not None:
                assert 0 < result[group]["mean"] < 1
            if result[group]["median"] is not None:
                assert 0 < result[group]["median"] < 1


class TestPlot:
    def test_plot_creates_file(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        with tempfile.NamedTemporaryFile(
            suffix=".png", delete=False
        ) as f:
            tmp_path = f.name
        try:
            args = Namespace(
                tree=TREE_SIMPLE, gene_trees=GENE_TREES,
                verbose=False, json=False, plot_output=tmp_path,
            )
            svc = EvoTempoMap(args)
            svc.run()
            assert os.path.exists(tmp_path)
            assert os.path.getsize(tmp_path) > 0
        finally:
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)

    def test_plot_not_created_without_flag(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=False, json=False, plot_output=None,
        )
        svc = EvoTempoMap(args)
        svc.run()  # Should not crash even though no plot output


class TestRun:
    def test_run_text_output(self, capsys):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap
        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=False, json=False, plot_output=None,
        )
        svc = EvoTempoMap(args)
        svc.run()
        captured = capsys.readouterr()
        assert "n_conc" in captured.out
        assert "n_disc" in captured.out
        assert "Global treeness" in captured.out
        assert "Branches tested" in captured.out

    def test_run_json_output(self, capsys):
        import json
        from phykit.services.tree.evo_tempo_map import EvoTempoMap
        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=False, json=True, plot_output=None,
        )
        svc = EvoTempoMap(args)
        svc.run()
        captured = capsys.readouterr()
        data = json.loads(captured.out)
        assert "branches" in data
        assert "global_" in data
        assert "n_gene_trees" in data["global_"]
        for branch in data["branches"]:
            assert "split" in branch
            assert "n_concordant" in branch
            assert "mann_whitney_p" in branch
            assert "fdr_p" in branch

    def test_run_verbose_output(self, capsys):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap
        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=True, json=False, plot_output=None,
        )
        svc = EvoTempoMap(args)
        svc.run()
        captured = capsys.readouterr()
        # Verbose should print per-branch raw lengths
        assert "Concordant lengths" in captured.out or "concordant" in captured.out.lower()


class TestCLIIntegration:
    def test_cli_runs_without_error(self, capsys):
        from phykit.phykit import Phykit

        testargs = [
            "phykit",
            "evo_tempo_map",
            "-t", TREE_SIMPLE,
            "-g", GENE_TREES,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        captured = capsys.readouterr()
        assert "n_conc" in captured.out

    def test_cli_alias_etm(self, capsys):
        from phykit.phykit import Phykit

        testargs = [
            "phykit",
            "etm",
            "-t", TREE_SIMPLE,
            "-g", GENE_TREES,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        captured = capsys.readouterr()
        assert "n_conc" in captured.out

    def test_cli_json_flag(self, capsys):
        import json
        from phykit.phykit import Phykit

        testargs = [
            "phykit",
            "evo_tempo_map",
            "-t", TREE_SIMPLE,
            "-g", GENE_TREES,
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        captured = capsys.readouterr()
        data = json.loads(captured.out)
        assert "branches" in data
