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
        assert len(parts) == 2
        float(parts[0])  # K
        float(parts[1])  # p_value

    def test_lambda_text_output(self, lambda_args, capsys):
        svc = PhylogeneticSignal(lambda_args)
        svc.run()
        out, _ = capsys.readouterr()
        parts = out.strip().split("\t")
        assert len(parts) == 3
        float(parts[0])  # lambda
        float(parts[1])  # log_likelihood
        float(parts[2])  # p_value

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
