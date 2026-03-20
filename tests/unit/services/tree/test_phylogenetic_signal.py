import pytest
import json
import numpy as np
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch

from phykit.services.tree.phylogenetic_signal import PhylogeneticSignal
import phykit.services.tree.phylogenetic_signal as ps_module
from phykit.errors import PhykitUserError


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_traits.tsv")
MULTI_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_multi_traits.tsv")


@pytest.fixture
def default_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=TRAITS_FILE,
        method="blombergs_k",
        permutations=1000,
        json=False,
    )


@pytest.fixture
def lambda_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=TRAITS_FILE,
        method="lambda",
        permutations=1000,
        json=False,
    )


class TestProcessArgs:
    def test_defaults(self):
        args = Namespace(tree="t.tre", trait_data="d.tsv")
        svc = PhylogeneticSignal(args)
        assert svc.tree_file_path == "t.tre"
        assert svc.trait_data_path == "d.tsv"
        assert svc.method == "blombergs_k"
        assert svc.permutations == 1000
        assert svc.json_output is False

    def test_json_true(self):
        args = Namespace(tree="t.tre", trait_data="d.tsv", json=True, method="lambda", permutations=500)
        svc = PhylogeneticSignal(args)
        assert svc.json_output is True
        assert svc.method == "lambda"
        assert svc.permutations == 500


class TestParseTraitFile:
    def test_valid_file(self, default_args):
        svc = PhylogeneticSignal(default_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        traits = svc._parse_trait_file(TRAITS_FILE, tips)
        assert len(traits) == 8
        assert traits["raccoon"] == pytest.approx(1.04)
        assert traits["weasel"] == pytest.approx(-0.30)

    def test_missing_file(self, default_args):
        svc = PhylogeneticSignal(default_args)
        with pytest.raises(PhykitUserError):
            svc._parse_trait_file("/nonexistent/path.tsv", ["a", "b", "c"])

    def test_non_numeric_value(self, default_args, tmp_path):
        trait_file = tmp_path / "bad.tsv"
        trait_file.write_text("taxon1\tabc\n")
        svc = PhylogeneticSignal(default_args)
        with pytest.raises(PhykitUserError):
            svc._parse_trait_file(str(trait_file), ["taxon1"])

    def test_wrong_columns(self, default_args, tmp_path):
        trait_file = tmp_path / "bad.tsv"
        trait_file.write_text("taxon1\t1.0\textra\n")
        svc = PhylogeneticSignal(default_args)
        with pytest.raises(PhykitUserError):
            svc._parse_trait_file(str(trait_file), ["taxon1"])

    def test_comments_and_blanks(self, default_args, tmp_path):
        trait_file = tmp_path / "good.tsv"
        trait_file.write_text("# comment\n\ntaxon1\t1.0\ntaxon2\t2.0\ntaxon3\t3.0\n")
        svc = PhylogeneticSignal(default_args)
        traits = svc._parse_trait_file(str(trait_file), ["taxon1", "taxon2", "taxon3"])
        assert len(traits) == 3

    def test_taxon_mismatch_warns(self, default_args, tmp_path, capsys):
        trait_file = tmp_path / "partial.tsv"
        trait_file.write_text("taxon1\t1.0\ntaxon2\t2.0\ntaxon3\t3.0\nextra\t4.0\n")
        svc = PhylogeneticSignal(default_args)
        traits = svc._parse_trait_file(
            str(trait_file), ["taxon1", "taxon2", "taxon3", "taxon4"]
        )
        _, err = capsys.readouterr()
        assert "Warning" in err
        assert len(traits) == 3

    def test_too_few_shared_taxa(self, default_args, tmp_path):
        trait_file = tmp_path / "few.tsv"
        trait_file.write_text("taxon1\t1.0\ntaxon2\t2.0\n")
        svc = PhylogeneticSignal(default_args)
        with pytest.raises(PhykitUserError):
            svc._parse_trait_file(str(trait_file), ["taxon1", "taxon2"])


class TestValidateTree:
    def test_too_few_tips(self, default_args, mocker):
        svc = PhylogeneticSignal(default_args)
        from Bio import Phylo
        from io import StringIO
        tree = Phylo.read(StringIO("(A:1.0,B:2.0);"), "newick")
        with pytest.raises(PhykitUserError):
            svc._validate_tree(tree)

    def test_missing_branch_lengths(self, default_args):
        svc = PhylogeneticSignal(default_args)
        from Bio import Phylo
        from io import StringIO
        tree = Phylo.read(StringIO("(A,B,C);"), "newick")
        with pytest.raises(PhykitUserError):
            svc._validate_tree(tree)

    def test_valid_tree_passes(self, default_args):
        svc = PhylogeneticSignal(default_args)
        tree = svc.read_tree_file()
        svc._validate_tree(tree)


class TestBuildVCVMatrix:
    def test_symmetric(self, default_args):
        svc = PhylogeneticSignal(default_args)
        tree = svc.read_tree_file()
        tips = sorted(svc.get_tip_names_from_tree(tree))
        vcv = svc._build_vcv_matrix(tree, tips)
        np.testing.assert_array_almost_equal(vcv, vcv.T)

    def test_diagonal_is_root_to_tip(self, default_args):
        svc = PhylogeneticSignal(default_args)
        tree = svc.read_tree_file()
        tips = sorted(svc.get_tip_names_from_tree(tree))
        vcv = svc._build_vcv_matrix(tree, tips)
        for i, name in enumerate(tips):
            expected = tree.distance(tree.root, name)
            assert vcv[i, i] == pytest.approx(expected, rel=1e-6)

    def test_correct_shape(self, default_args):
        svc = PhylogeneticSignal(default_args)
        tree = svc.read_tree_file()
        tips = sorted(svc.get_tip_names_from_tree(tree))
        vcv = svc._build_vcv_matrix(tree, tips)
        assert vcv.shape == (8, 8)


class TestBlombergsK:
    def test_returns_expected_keys(self, default_args):
        svc = PhylogeneticSignal(default_args)
        tree = svc.read_tree_file()
        tips = sorted(svc.get_tip_names_from_tree(tree))
        traits = svc._parse_trait_file(TRAITS_FILE, tips)
        ordered = sorted(traits.keys())
        x = np.array([traits[n] for n in ordered])
        vcv = svc._build_vcv_matrix(tree, ordered)
        result = svc._blombergs_k(x, vcv, 100)
        assert "K" in result
        assert "p_value" in result
        assert "permutations" in result
        assert result["permutations"] == 100

    def test_k_matches_phytools(self, default_args):
        """K must match R phytools::phylosig(method='K') reference value."""
        svc = PhylogeneticSignal(default_args)
        tree = svc.read_tree_file()
        tips = sorted(svc.get_tip_names_from_tree(tree))
        traits = svc._parse_trait_file(TRAITS_FILE, tips)
        ordered = sorted(traits.keys())
        x = np.array([traits[n] for n in ordered])
        vcv = svc._build_vcv_matrix(tree, ordered)
        result = svc._blombergs_k(x, vcv, 1000)
        # R reference: phylosig(tree, x, method="K") = 0.584216
        assert result["K"] == pytest.approx(0.584216, abs=1e-4)
        assert 0 <= result["p_value"] <= 1


class TestPagelsLambda:
    def test_returns_expected_keys(self, default_args):
        svc = PhylogeneticSignal(default_args)
        tree = svc.read_tree_file()
        tips = sorted(svc.get_tip_names_from_tree(tree))
        traits = svc._parse_trait_file(TRAITS_FILE, tips)
        ordered = sorted(traits.keys())
        x = np.array([traits[n] for n in ordered])
        vcv = svc._build_vcv_matrix(tree, ordered)
        result = svc._pagels_lambda(x, vcv)
        assert "lambda" in result
        assert "log_likelihood" in result
        assert "p_value" in result

    def test_lambda_matches_phytools(self, default_args):
        """Lambda and logL must match R phytools::phylosig(method='lambda') reference."""
        svc = PhylogeneticSignal(default_args)
        tree = svc.read_tree_file()
        tips = sorted(svc.get_tip_names_from_tree(tree))
        traits = svc._parse_trait_file(TRAITS_FILE, tips)
        ordered = sorted(traits.keys())
        x = np.array([traits[n] for n in ordered])
        vcv = svc._build_vcv_matrix(tree, ordered)
        result = svc._pagels_lambda(x, vcv)
        # R reference: phylosig(tree, x, method="lambda", test=TRUE)
        #   lambda = 0.999927, logL = -11.5697, P = 0.716569
        assert result["lambda"] == pytest.approx(0.9999, abs=1e-3)
        assert result["log_likelihood"] == pytest.approx(-11.5697, abs=1e-3)
        assert result["p_value"] == pytest.approx(0.7166, abs=1e-2)
        assert 0 <= result["lambda"] <= 1


class TestRun:
    def test_blombergs_k_text_output(self, default_args, capsys):
        svc = PhylogeneticSignal(default_args)
        svc.run()
        out, _ = capsys.readouterr()
        parts = out.strip().split("\t")
        assert len(parts) == 3
        float(parts[0])  # K
        float(parts[1])  # p_value
        float(parts[2])  # r_squared_phylo

    def test_lambda_text_output(self, lambda_args, capsys):
        svc = PhylogeneticSignal(lambda_args)
        svc.run()
        out, _ = capsys.readouterr()
        parts = out.strip().split("\t")
        assert len(parts) == 4
        float(parts[0])  # lambda
        float(parts[1])  # log_likelihood
        float(parts[2])  # p_value
        float(parts[3])  # r_squared_phylo

    def test_json_output_blombergs_k(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            method="blombergs_k",
            permutations=100,
            json=True,
        )
        svc = PhylogeneticSignal(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "K" in data
        assert "p_value" in data
        assert "permutations" in data

    def test_json_output_lambda(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            method="lambda",
            permutations=1000,
            json=True,
        )
        svc = PhylogeneticSignal(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "lambda" in data
        assert "log_likelihood" in data
        assert "p_value" in data


GENE_TREES_FILE = str(SAMPLE_FILES / "gene_trees_simple.nwk")


class TestDiscordanceVCV:
    def test_run_with_gene_trees(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            method="blombergs_k",
            permutations=100,
            json=True,
            gene_trees=GENE_TREES_FILE,
        )
        svc = PhylogeneticSignal(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "K" in data
        assert "p_value" in data
        assert "vcv_metadata" in data
        assert data["vcv_metadata"]["n_gene_trees"] == 10

    def test_all_concordant_gene_trees(self, capsys):
        """When all gene trees are identical to species tree, K should be
        very close to species-tree-only K."""
        # Run without gene trees
        args_no_gt = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            method="blombergs_k",
            permutations=100,
            json=True,
        )
        svc = PhylogeneticSignal(args_no_gt)
        svc.run()
        out, _ = capsys.readouterr()
        k_no_gt = json.loads(out)["K"]

        # Run with species tree as gene trees (all concordant)
        import tempfile, shutil
        from Bio import Phylo
        tree = Phylo.read(TREE_SIMPLE, "newick")
        with tempfile.NamedTemporaryFile(mode="w", suffix=".nwk", delete=False) as f:
            for _ in range(3):
                Phylo.write(tree, f, "newick")
            concordant_path = f.name

        try:
            args_gt = Namespace(
                tree=TREE_SIMPLE,
                trait_data=TRAITS_FILE,
                method="blombergs_k",
                permutations=100,
                json=True,
                gene_trees=concordant_path,
            )
            svc = PhylogeneticSignal(args_gt)
            svc.run()
            out, _ = capsys.readouterr()
            k_gt = json.loads(out)["K"]
            assert k_gt == pytest.approx(k_no_gt, rel=0.01)
        finally:
            import os
            os.unlink(concordant_path)

    def test_discordance_vcv_differs_from_species(self, capsys):
        """With discordant gene trees, K should differ from species-tree-only K."""
        args_no_gt = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            method="blombergs_k",
            permutations=100,
            json=True,
        )
        svc = PhylogeneticSignal(args_no_gt)
        svc.run()
        out, _ = capsys.readouterr()
        k_no_gt = json.loads(out)["K"]

        args_gt = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            method="blombergs_k",
            permutations=100,
            json=True,
            gene_trees=GENE_TREES_FILE,
        )
        svc = PhylogeneticSignal(args_gt)
        svc.run()
        out, _ = capsys.readouterr()
        k_gt = json.loads(out)["K"]
        # They should differ (gene trees 7-9 are discordant)
        assert k_gt != pytest.approx(k_no_gt, abs=1e-6)

    def test_lambda_with_gene_trees(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            method="lambda",
            permutations=1000,
            json=True,
            gene_trees=GENE_TREES_FILE,
        )
        svc = PhylogeneticSignal(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "lambda" in data
        assert "vcv_metadata" in data


class TestEffectSize:
    def test_r2_phylo_blombergs_k(self, capsys):
        """JSON output for blombergs_k contains r_squared_phylo as a finite float."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            method="blombergs_k",
            permutations=100,
            json=True,
        )
        svc = PhylogeneticSignal(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "r_squared_phylo" in data
        assert isinstance(data["r_squared_phylo"], float)
        assert np.isfinite(data["r_squared_phylo"])

    def test_r2_phylo_lambda(self, capsys):
        """JSON output for lambda contains r_squared_phylo as a finite float."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            method="lambda",
            permutations=1000,
            json=True,
        )
        svc = PhylogeneticSignal(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "r_squared_phylo" in data
        assert isinstance(data["r_squared_phylo"], float)
        assert np.isfinite(data["r_squared_phylo"])

    def test_r2_phylo_positive_for_phylogenetic_trait(self, capsys):
        """Body mass trait has strong phylogenetic signal, so R2_phylo > 0."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            method="blombergs_k",
            permutations=100,
            json=True,
        )
        svc = PhylogeneticSignal(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert data["r_squared_phylo"] > 0

    def test_r2_phylo_in_text_output(self, capsys):
        """Text output includes R2_phylo value."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            method="blombergs_k",
            permutations=100,
            json=False,
        )
        svc = PhylogeneticSignal(args)
        svc.run()
        out, _ = capsys.readouterr()
        parts = out.strip().split("\t")
        # Third column is r_squared_phylo
        assert len(parts) == 3
        r2 = float(parts[2])
        assert np.isfinite(r2)

    def test_r2_phylo_matches_r_reference(self, capsys):
        """R²_phylo must match R phytools/manual GLS reference value.

        R reference (tests/r_validation/validate_signal_r2.R):
          sigma2_bm_manual  = 0.0384065703
          sigma2_wn_manual  = 0.7667234375
          r2_phylo_manual   = 0.9499081827
        """
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            method="blombergs_k",
            permutations=100,
            json=True,
        )
        svc = PhylogeneticSignal(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert data["r_squared_phylo"] == pytest.approx(0.9499081827, abs=1e-6)

    def test_r2_phylo_negative_for_anti_phylogenetic_trait(self, tmp_path, capsys):
        """R²_phylo < 0 when closely related species have very different values.

        The design doc states: 'R²_phylo < 0: possible (phylogeny actively
        misleads), valid. Report as-is.' This test verifies that negative R²
        values are correctly computed and not clipped to 0 or NaN.

        We use a balanced tree with long internal branches (strong expected
        phylogenetic signal) but assign opposite trait values to sister taxa,
        maximally violating the BM assumption.
        """
        # Tree with strong phylogenetic structure (long shared paths)
        tree_file = tmp_path / "balanced.tre"
        tree_file.write_text("((A:0.1,B:0.1):10,(C:0.1,D:0.1):10);")

        # Anti-phylogenetic: sisters have opposite values
        trait_file = tmp_path / "anti_phylo_traits.tsv"
        trait_file.write_text("A\t10.0\nB\t-10.0\nC\t5.0\nD\t-5.0\n")

        args = Namespace(
            tree=str(tree_file),
            trait_data=str(trait_file),
            method="blombergs_k",
            permutations=100,
            json=True,
        )
        svc = PhylogeneticSignal(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        # R²_phylo should be strongly negative: phylogeny misleads
        assert data["r_squared_phylo"] < 0, (
            f"Expected negative R²_phylo for anti-phylogenetic trait, "
            f"got {data['r_squared_phylo']}"
        )
        # Should be a real number, not NaN
        assert np.isfinite(data["r_squared_phylo"])
        # Reference: R²_phylo = -9.0 for this exact setup
        assert data["r_squared_phylo"] == pytest.approx(-9.0, abs=0.1)


class TestKmult:
    def test_kmult_returns_positive(self):
        """K_mult should be > 0 for real trait data."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="blombergs_k",
            permutations=100,
            json=True,
            multivariate=True,
        )
        svc = PhylogeneticSignal(args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits_multi = svc._parse_multi_trait_file(MULTI_TRAITS_FILE, tips)
        ordered = sorted(traits_multi.keys())
        from phykit.services.tree.vcv_utils import build_vcv_matrix
        vcv = build_vcv_matrix(tree, ordered)
        p = len(trait_names)
        Y = np.array([[traits_multi[name][j] for j in range(p)] for name in ordered])
        result = svc._kmult(Y, vcv, 100)
        assert result["K_mult"] > 0

    def test_kmult_permutation_pvalue(self):
        """p-value from K_mult permutation test must be between 0 and 1."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="blombergs_k",
            permutations=100,
            json=True,
            multivariate=True,
        )
        svc = PhylogeneticSignal(args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits_multi = svc._parse_multi_trait_file(MULTI_TRAITS_FILE, tips)
        ordered = sorted(traits_multi.keys())
        from phykit.services.tree.vcv_utils import build_vcv_matrix
        vcv = build_vcv_matrix(tree, ordered)
        p = len(trait_names)
        Y = np.array([[traits_multi[name][j] for j in range(p)] for name in ordered])
        result = svc._kmult(Y, vcv, 100)
        assert 0 <= result["p_value"] <= 1

    def test_kmult_single_trait_matches_k(self):
        """K_mult with 1 trait should approximately equal univariate K."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            method="blombergs_k",
            permutations=100,
            json=True,
            multivariate=False,
        )
        svc = PhylogeneticSignal(args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        traits = svc._parse_trait_file(TRAITS_FILE, tips)
        ordered = sorted(traits.keys())
        x = np.array([traits[n] for n in ordered])
        from phykit.services.tree.vcv_utils import build_vcv_matrix
        vcv = build_vcv_matrix(tree, ordered)

        # Univariate K
        result_k = svc._blombergs_k(x, vcv, 100)

        # K_mult with single trait (n x 1 matrix)
        Y = x.reshape(-1, 1)
        result_kmult = svc._kmult(Y, vcv, 100)

        assert result_kmult["K_mult"] == pytest.approx(result_k["K"], rel=1e-6)

    def test_multivariate_flag_uses_multi_trait_file(self, capsys):
        """--multivariate flag reads multi-column TSV and runs K_mult."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="blombergs_k",
            permutations=100,
            json=True,
            multivariate=True,
        )
        svc = PhylogeneticSignal(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "K_mult" in data
        assert "p_value" in data
        assert "n_traits" in data
        assert data["n_traits"] == 3
        assert "permutations" in data
        assert data["permutations"] == 100

    def test_multivariate_text_output(self, capsys):
        """Text output for K_mult should have 4 tab-separated fields."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="blombergs_k",
            permutations=100,
            json=False,
            multivariate=True,
        )
        svc = PhylogeneticSignal(args)
        svc.run()
        out, _ = capsys.readouterr()
        parts = out.strip().split("\t")
        assert len(parts) == 4
        float(parts[0])  # K_mult
        float(parts[1])  # p_value
        int(parts[2])    # n_traits
        int(parts[3])    # permutations

    def test_multivariate_with_lambda_raises_error(self):
        """Using --multivariate with --method lambda should raise an error."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="lambda",
            permutations=100,
            json=False,
            multivariate=True,
        )
        svc = PhylogeneticSignal(args)
        with pytest.raises(PhykitUserError):
            svc.run()

    def test_multivariate_json_output(self, capsys):
        """JSON output for K_mult should contain all expected fields."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="blombergs_k",
            permutations=100,
            json=True,
            multivariate=True,
        )
        svc = PhylogeneticSignal(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "K_mult" in data
        assert isinstance(data["K_mult"], float)
        assert "p_value" in data
        assert isinstance(data["p_value"], float)
        assert "n_traits" in data
        assert data["n_traits"] == 3
        assert "permutations" in data
        assert data["permutations"] == 100

    def test_multivariate_with_gene_trees(self, capsys):
        """K_mult should work with discordance-aware VCV (--gene-trees)."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="blombergs_k",
            permutations=100,
            json=True,
            multivariate=True,
            gene_trees=str(SAMPLE_FILES / "gene_trees_simple.nwk"),
        )
        svc = PhylogeneticSignal(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "K_mult" in data
        assert "vcv_metadata" in data
