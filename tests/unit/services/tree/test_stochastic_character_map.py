import json
import os
import subprocess
import sys
import tempfile
from io import StringIO

import pytest
import numpy as np
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, TreeMixin

import phykit.services.tree.stochastic_character_map as scm_module
from phykit.services.tree.stochastic_character_map import StochasticCharacterMap
from phykit.errors import PhykitUserError


def test_module_import_does_not_import_numpy():
    code = """
import sys
import phykit.services.tree.stochastic_character_map as module
assert hasattr(module.np, "__getattr__")
assert callable(module.print_json)
assert "typing" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
assert "phykit.helpers.discrete_models" not in sys.modules
assert "phykit.helpers.circular_layout" not in sys.modules
assert "phykit.helpers.color_annotations" not in sys.modules
assert "numpy" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_lazy_numpy_caches_resolved_attributes():
    lazy_np = scm_module._LazyNumpy()

    zeros_attr = lazy_np.zeros

    assert lazy_np.__dict__["zeros"] is zeros_attr
    assert lazy_np.zeros is zeros_attr
    assert lazy_np._module is not None


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
DISCRETE_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_discrete_traits.tsv")


@pytest.fixture
def default_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=DISCRETE_TRAITS_FILE,
        trait="diet",
        model="ER",
        nsim=10,
        seed=42,
        plot=None,
        json=False,
    )


@pytest.fixture
def ard_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=DISCRETE_TRAITS_FILE,
        trait="diet",
        model="ARD",
        nsim=10,
        seed=42,
        plot=None,
        json=False,
    )


@pytest.fixture
def json_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=DISCRETE_TRAITS_FILE,
        trait="diet",
        model="ER",
        nsim=10,
        seed=42,
        plot=None,
        json=True,
    )


class TestProcessArgs:
    def test_defaults(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            trait="diet",
        )
        svc = StochasticCharacterMap.__new__(StochasticCharacterMap)
        parsed = svc.process_args(args)
        assert parsed["tree_file_path"] == "t.tre"
        assert parsed["trait_data_path"] == "d.tsv"
        assert parsed["trait_column"] == "diet"
        assert parsed["model"] == "ER"
        assert parsed["nsim"] == 100
        assert parsed["seed"] is None
        assert parsed["plot_output"] is None
        assert parsed["json_output"] is False

    def test_overrides(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            trait="diet",
            model="ARD",
            nsim=500,
            seed=123,
            plot="out.png",
            json=True,
        )
        svc = StochasticCharacterMap.__new__(StochasticCharacterMap)
        parsed = svc.process_args(args)
        assert parsed["model"] == "ARD"
        assert parsed["nsim"] == 500
        assert parsed["seed"] == 123
        assert parsed["plot_output"] == "out.png"
        assert parsed["json_output"] is True


class TestDiscreteTraitParsing:
    def test_basic_parsing(self, default_args):
        svc = StochasticCharacterMap(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        tip_states = svc._parse_discrete_trait_file(
            DISCRETE_TRAITS_FILE, "diet", tree_tips
        )
        assert len(tip_states) == 8
        assert tip_states["raccoon"] == "carnivore"
        assert tip_states["monkey"] == "herbivore"
        assert tip_states["dog"] == "omnivore"

    def test_missing_column(self, default_args):
        svc = StochasticCharacterMap(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        with pytest.raises(PhykitUserError):
            svc._parse_discrete_trait_file(
                DISCRETE_TRAITS_FILE, "nonexistent", tree_tips
            )

    def test_comments_and_blanks(self, default_args):
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f:
            f.write("# This is a comment\n")
            f.write("\n")
            f.write("taxon\ttrait\n")
            f.write("raccoon\tA\n")
            f.write("bear\tA\n")
            f.write("# Another comment\n")
            f.write("sea_lion\tB\n")
            f.write("seal\tB\n")
            f.write("monkey\tA\n")
            f.write("cat\tB\n")
            f.write("weasel\tA\n")
            f.write("dog\tB\n")
            tmppath = f.name

        try:
            svc = StochasticCharacterMap(default_args)
            tree = svc.read_tree_file()
            tree_tips = svc.get_tip_names_from_tree(tree)
            tip_states = svc._parse_discrete_trait_file(tmppath, "trait", tree_tips)
            assert len(tip_states) == 8
        finally:
            os.unlink(tmppath)


class TestTreePruning:
    def test_count_missing_tip_states_uses_direct_traversal(self, monkeypatch):
        tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("generic terminal traversal should not be used")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        missing_count = StochasticCharacterMap._count_missing_tip_states(
            tree,
            {"A": "x", "C": "y"},
        )

        assert missing_count == 2

    def test_count_missing_tip_states_handles_mixed_child_counts(self, monkeypatch):
        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);"),
            "newick",
        )

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("standard tree counter should use direct traversal")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        missing_count = StochasticCharacterMap._count_missing_tip_states(
            tree,
            {"A": "x", "C": "y", "F": "z"},
        )

        assert missing_count == 3

    def test_prune_tree_to_tip_states_uses_batch_prune_for_standard_tree(self, monkeypatch):
        from Bio import Phylo
        from Bio.Phylo.BaseTree import TreeMixin
        from io import StringIO

        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1,G:1);"),
            "newick",
        )

        def fail_prune(*_args, **_kwargs):
            raise AssertionError("per-tip prune should not be used")

        monkeypatch.setattr(TreeMixin, "prune", fail_prune)

        removed = StochasticCharacterMap._prune_tree_to_tip_states(
            tree,
            {"A": "x", "B": "x", "D": "y", "F": "y"},
        )

        assert removed == 3
        assert {tip.name for tip in tree.get_terminals()} == {"A", "B", "D", "F"}

    def test_prune_tree_to_tip_states_prunes_terminal_objects(self, monkeypatch):
        from Bio import Phylo
        from Bio.Phylo.BaseTree import TreeMixin
        from io import StringIO

        tree = Phylo.read(StringIO("(A:1,B:1);"), "newick")
        original_prune = TreeMixin.prune

        def assert_object_prune(self, target=None, **kwargs):
            assert not isinstance(target, str)
            return original_prune(self, target=target, **kwargs)

        monkeypatch.setattr(TreeMixin, "prune", assert_object_prune)

        removed = StochasticCharacterMap._prune_tree_to_tip_states(
            tree,
            {"A": "x"},
        )

        assert removed == 1
        assert {tip.name for tip in tree.get_terminals()} == {"A"}


class TestTraversalHelpers:
    def test_build_parent_map_uses_direct_traversal(self, default_args, monkeypatch):
        tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )
        svc = StochasticCharacterMap(default_args)

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("generic traversal should not be used")

        monkeypatch.setattr(tree, "find_clades", fail_find_clades)

        parent_map = svc._build_parent_map(tree)

        assert parent_map[id(tree.root.clades[0])] is tree.root
        assert parent_map[id(tree.root.clades[1])] is tree.root
        assert parent_map[id(tree.root.clades[0].clades[0])] is tree.root.clades[0]
        assert parent_map[id(tree.root.clades[1].clades[1])] is tree.root.clades[1]

    def test_build_parent_map_handles_mixed_child_counts(
        self, default_args, monkeypatch
    ):
        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);"),
            "newick",
        )
        svc = StochasticCharacterMap(default_args)

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("generic traversal should not be used")

        monkeypatch.setattr(tree, "find_clades", fail_find_clades)

        parent_map = svc._build_parent_map(tree)
        binary = tree.root.clades[1]
        trifurcation = tree.root.clades[2]

        assert parent_map[id(tree.root.clades[0])] is tree.root
        assert parent_map[id(binary)] is tree.root
        assert parent_map[id(trifurcation)] is tree.root
        assert parent_map[id(binary.clades[0])] is binary
        assert parent_map[id(binary.clades[1])] is binary
        assert parent_map[id(trifurcation.clades[0])] is trifurcation
        assert parent_map[id(trifurcation.clades[1])] is trifurcation
        assert parent_map[id(trifurcation.clades[2])] is trifurcation

    def test_build_simulation_metadata_uses_direct_preorder_traversal(
        self, default_args, monkeypatch
    ):
        tree = Phylo.read(
            StringIO("((A:1,B:2):3,(C:4,D:5):6);"),
            "newick",
        )
        svc = StochasticCharacterMap(default_args)
        tip_states = {"A": "x", "B": "y", "C": "x", "D": "y"}
        states = ["x", "y"]
        parent_map = svc._build_parent_map(tree)

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("generic traversal should not be used")

        monkeypatch.setattr(tree, "find_clades", fail_find_clades)

        simulation_nodes, branch_nodes = svc._build_simulation_metadata(
            tree,
            tip_states,
            states,
            parent_map,
        )

        left, right = tree.root.clades
        expected_simulation_nodes = [
            (id(left), id(tree.root), 3.0, None),
            (id(left.clades[0]), id(left), 1.0, 0),
            (id(left.clades[1]), id(left), 2.0, 1),
            (id(right), id(tree.root), 6.0, None),
            (id(right.clades[0]), id(right), 4.0, 0),
            (id(right.clades[1]), id(right), 5.0, 1),
        ]
        expected_branch_nodes = [
            (clade_id, parent_id, branch_length)
            for clade_id, parent_id, branch_length, _ in expected_simulation_nodes
        ]

        assert simulation_nodes == expected_simulation_nodes
        assert branch_nodes == expected_branch_nodes

    def test_build_simulation_metadata_reuses_optimized_preorder(
        self, default_args
    ):
        class NoReversedList(list):
            def __reversed__(self):
                raise AssertionError("metadata traversal should not call reversed")

        svc = StochasticCharacterMap(default_args)
        root = Clade(name="root")
        left = Clade(branch_length=3.0, name="left")
        right = Clade(branch_length=4.0, name="right")
        left.clades = NoReversedList(
            [
                Clade(branch_length=1.0, name="A"),
                Clade(branch_length=2.0, name="B"),
            ]
        )
        root.clades = NoReversedList([left, right])
        tree = type("Tree", (), {"root": root})()
        parent_map = {
            id(left): root,
            id(right): root,
            id(left.clades[0]): left,
            id(left.clades[1]): left,
        }

        simulation_nodes, branch_nodes = svc._build_simulation_metadata(
            tree,
            {"A": "x", "B": "y", "right": "x"},
            ["x", "y"],
            parent_map,
        )

        expected_simulation_nodes = [
            (id(left), id(root), 3.0, None),
            (id(left.clades[0]), id(left), 1.0, 0),
            (id(left.clades[1]), id(left), 2.0, 1),
            (id(right), id(root), 4.0, 0),
        ]
        assert simulation_nodes == expected_simulation_nodes
        assert branch_nodes == [
            (clade_id, parent_id, branch_length)
            for clade_id, parent_id, branch_length, _ in expected_simulation_nodes
        ]


class TestQMatrixFitting:
    """Validated against R 4.4.0 with phytools::fitMk on same tree and traits."""

    def test_er_model_loglik(self, default_args):
        svc = StochasticCharacterMap(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        tip_states = svc._parse_discrete_trait_file(
            DISCRETE_TRAITS_FILE, "diet", tree_tips
        )
        states = sorted(set(tip_states.values()))
        Q, loglik = svc._fit_q_matrix(tree, tip_states, states, "ER")

        # R reference: fitMk(tree, trait, model="ER") -> loglik = -8.788898
        # Python finds a marginally better point on a flat likelihood surface
        # (-8.787360 vs -8.788898), both valid ML estimates
        assert loglik == pytest.approx(-8.7889, abs=0.01)

        # Q should have equal off-diagonal rates for ER
        k = len(states)
        off_diag = []
        for i in range(k):
            for j in range(k):
                if i != j:
                    off_diag.append(Q[i, j])
        # All off-diagonal rates should be equal for ER
        assert all(abs(r - off_diag[0]) < 1e-6 for r in off_diag)
        # Rate should be positive
        assert off_diag[0] > 0

    def test_sym_model(self, default_args):
        svc = StochasticCharacterMap(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        tip_states = svc._parse_discrete_trait_file(
            DISCRETE_TRAITS_FILE, "diet", tree_tips
        )
        states = sorted(set(tip_states.values()))

        Q_er, loglik_er = svc._fit_q_matrix(tree, tip_states, states, "ER")
        Q_sym, loglik_sym = svc._fit_q_matrix(tree, tip_states, states, "SYM")

        # SYM should give same or better loglik than ER (superset model)
        assert loglik_sym >= loglik_er - 0.01

        # Q should be symmetric
        k = len(states)
        for i in range(k):
            for j in range(k):
                if i != j:
                    assert Q_sym[i, j] == pytest.approx(Q_sym[j, i], abs=1e-4)

    def test_ard_model_loglik(self, default_args):
        svc = StochasticCharacterMap(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        tip_states = svc._parse_discrete_trait_file(
            DISCRETE_TRAITS_FILE, "diet", tree_tips
        )
        states = sorted(set(tip_states.values()))

        _, loglik_er = svc._fit_q_matrix(tree, tip_states, states, "ER")
        Q_ard, loglik_ard = svc._fit_q_matrix(tree, tip_states, states, "ARD")

        # R reference: fitMk(tree, trait, model="ARD") -> loglik = -8.430467
        # Python's multi-start optimizer finds a better local optimum (~-8.385)
        # on this 6-parameter surface. Both are valid ML estimates.
        assert loglik_ard >= -8.45  # at least as good as R
        assert np.isfinite(loglik_ard)

        # ARD should give same or better loglik than ER (superset model)
        assert loglik_ard >= loglik_er - 0.01

    def test_q_matrix_structure(self, default_args):
        svc = StochasticCharacterMap(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        tip_states = svc._parse_discrete_trait_file(
            DISCRETE_TRAITS_FILE, "diet", tree_tips
        )
        states = sorted(set(tip_states.values()))
        Q, _ = svc._fit_q_matrix(tree, tip_states, states, "ER")

        k = len(states)
        # Rows should sum to 0
        for i in range(k):
            assert abs(Q[i, :].sum()) < 1e-10

        # Off-diagonal should be positive
        for i in range(k):
            for j in range(k):
                if i != j:
                    assert Q[i, j] > 0

        # Diagonal should be negative
        for i in range(k):
            assert Q[i, i] < 0


class TestFelsensteinPruning:
    def test_conditional_likelihoods(self, default_args):
        svc = StochasticCharacterMap(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        tip_states = svc._parse_discrete_trait_file(
            DISCRETE_TRAITS_FILE, "diet", tree_tips
        )
        states = sorted(set(tip_states.values()))
        k = len(states)
        pi = np.ones(k) / k

        # Use a simple Q
        Q = np.array([[-2, 1, 1], [1, -2, 1], [1, 1, -2]], dtype=float)
        cond_liks, loglik = svc._felsenstein_pruning(
            tree, tip_states, Q, pi, states
        )

        # Should have conditional likelihoods for all nodes
        node_count = sum(1 for _ in tree.find_clades())
        assert len(cond_liks) == node_count

        # Log-likelihood should be finite and negative
        assert np.isfinite(loglik)
        assert loglik < 0

        # Tip conditional likelihoods should be indicator vectors
        for tip in tree.get_terminals():
            lik = cond_liks[id(tip)]
            assert lik.sum() == pytest.approx(1.0)
            assert all(v in (0.0, 1.0) for v in lik)


class TestStochasticMapping:
    def test_sample_categorical_matches_numpy_choice_stream(self):
        probs = np.array([0.2, 0.3, 0.5])
        rng_choice = np.random.default_rng(123)
        rng_fast = np.random.default_rng(123)

        expected = [rng_choice.choice(len(probs), p=probs) for _ in range(100)]
        observed = [
            StochasticCharacterMap._sample_categorical(probs, rng_fast)
            for _ in range(100)
        ]

        assert observed == expected
        assert rng_fast.random() == pytest.approx(rng_choice.random())

    def test_sample_categorical_cdf_preserves_rng_stream(self):
        probs = np.array([0.2, 0.3, 0.5])
        cdf = np.cumsum(probs)
        rng_probs = np.random.default_rng(123)
        rng_cdf = np.random.default_rng(123)

        expected = [
            StochasticCharacterMap._sample_categorical(probs, rng_probs)
            for _ in range(100)
        ]
        observed = [
            StochasticCharacterMap._sample_categorical_cdf(cdf, rng_cdf)
            for _ in range(100)
        ]

        assert observed == expected
        assert rng_cdf.random() == pytest.approx(rng_probs.random())

    @pytest.mark.parametrize(
        "cdf,draws",
        [
            (np.array([0.2, 1.0]), [0.0, 0.2, 0.999, 1.0]),
            (np.array([0.2, 0.5, 1.0]), [0.0, 0.2, 0.499, 0.5, 0.999, 1.0]),
            (
                np.array([0.1, 0.3, 0.7, 1.0]),
                [0.0, 0.1, 0.299, 0.3, 0.699, 0.7, 0.999, 1.0],
            ),
            (
                np.array([0.1, 0.3, 0.5, 0.7, 1.0]),
                [0.0, 0.1, 0.299, 0.3, 0.499, 0.5, 0.999, 1.0],
            ),
            (
                np.array([0.05, 0.15, 0.35, 0.55, 0.75, 1.0]),
                [0.0, 0.05, 0.149, 0.15, 0.349, 0.35, 0.999, 1.0],
            ),
            (
                np.array([0.05, 0.1, 0.2, 0.35, 0.55, 0.8, 1.0]),
                [0.0, 0.05, 0.099, 0.1, 0.349, 0.35, 0.799, 0.8, 1.0],
            ),
            (
                np.array([0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 0.9, 1.0]),
                [0.0, 0.05, 0.099, 0.1, 0.499, 0.5, 0.899, 0.9, 0.999, 1.0],
            ),
        ],
    )
    def test_sample_categorical_cdf_small_state_boundaries(self, cdf, draws):
        class FixedRng:
            def __init__(self, value):
                self.value = value

            def random(self):
                return self.value

        for draw in draws:
            expected = int(np.searchsorted(cdf, draw, side="right"))
            if expected == len(cdf):
                expected -= 1

            observed = StochasticCharacterMap._sample_categorical_cdf(
                cdf,
                FixedRng(draw),
            )

            assert observed == expected

    def test_simulate_branch_history_reuses_precomputed_transition_probs(
        self, default_args
    ):
        svc = StochasticCharacterMap(default_args)
        Q = np.array(
            [
                [-0.7, 0.4, 0.3],
                [0.2, -0.5, 0.3],
                [0.1, 0.2, -0.3],
            ],
            dtype=float,
        )
        k = Q.shape[0]
        rates = -np.diag(Q)
        transition_probs = svc._transition_probs_by_state(Q, k)

        baseline = svc._simulate_branch_history(
            Q,
            0,
            2,
            1.5,
            k,
            np.random.default_rng(11),
        )
        optimized = svc._simulate_branch_history(
            Q,
            0,
            2,
            1.5,
            k,
            np.random.default_rng(11),
            rates=rates,
            transition_probs_by_state=transition_probs,
        )

        assert optimized == baseline
        assert np.allclose(transition_probs[0], [0.0, 4.0 / 7.0, 3.0 / 7.0])

    def test_simulate_branch_history_reuses_precomputed_transition_cdfs(
        self, default_args, monkeypatch
    ):
        svc = StochasticCharacterMap(default_args)
        Q = np.array(
            [
                [-0.7, 0.4, 0.3],
                [0.2, -0.5, 0.3],
                [0.1, 0.2, -0.3],
            ],
            dtype=float,
        )
        k = Q.shape[0]
        rates = -np.diag(Q)
        transition_probs = svc._transition_probs_by_state(Q, k)
        transition_cdfs = svc._transition_cdfs_by_state(transition_probs)

        baseline = svc._simulate_branch_history(
            Q,
            0,
            2,
            1.5,
            k,
            np.random.default_rng(11),
        )

        def fail_sample_categorical(*_args, **_kwargs):
            raise AssertionError("precomputed CDFs should be sampled directly")

        monkeypatch.setattr(svc, "_sample_categorical", fail_sample_categorical)

        optimized = svc._simulate_branch_history(
            Q,
            0,
            2,
            1.5,
            k,
            np.random.default_rng(11),
            rates=rates,
            transition_cdfs_by_state=transition_cdfs,
        )

        assert optimized == baseline

    def test_simulate_branch_history_uses_deterministic_two_state_targets(
        self, default_args, monkeypatch
    ):
        svc = StochasticCharacterMap(default_args)
        Q = np.array([[-0.7, 0.7], [0.4, -0.4]], dtype=float)
        k = Q.shape[0]
        rates, transition_cdfs, deterministic_targets = (
            svc._prepare_branch_history_context(Q, k)
        )

        baseline = svc._simulate_branch_history(
            Q,
            0,
            1,
            1.5,
            k,
            np.random.default_rng(11),
            rates=rates,
            transition_cdfs_by_state=transition_cdfs,
        )

        def fail_sample_categorical_cdf(*_args, **_kwargs):
            raise AssertionError("deterministic two-state transitions need no CDF search")

        monkeypatch.setattr(
            svc,
            "_sample_categorical_cdf",
            fail_sample_categorical_cdf,
        )

        optimized = svc._simulate_branch_history(
            Q,
            0,
            1,
            1.5,
            k,
            np.random.default_rng(11),
            rates=rates,
            transition_cdfs_by_state=transition_cdfs,
            deterministic_next_states=deterministic_targets,
        )

        assert deterministic_targets == [1, 0]
        assert optimized == baseline

    def test_prepare_branch_history_context_binary_fast_path(
        self, default_args, monkeypatch
    ):
        svc = StochasticCharacterMap(default_args)
        Q = np.array([[0.0, 0.7], [0.0, 0.0]], dtype=float)

        def fail_transition_probs(*_args, **_kwargs):
            raise AssertionError("binary context should not build transition probs")

        monkeypatch.setattr(
            svc,
            "_transition_probs_by_state",
            fail_transition_probs,
        )

        rates, transition_cdfs, deterministic_targets = (
            svc._prepare_branch_history_context(Q, Q.shape[0])
        )

        np.testing.assert_allclose(rates, -Q.diagonal())
        np.testing.assert_allclose(transition_cdfs[0], [0.0, 1.0])
        assert transition_cdfs[1] is None
        assert deterministic_targets == [1, -1]

    def test_prepare_branch_history_context_uses_diagonal_rate_view(
        self, default_args, monkeypatch
    ):
        svc = StochasticCharacterMap(default_args)
        Q = np.array(
            [
                [-0.7, 0.4, 0.3],
                [0.2, -0.5, 0.3],
                [0.1, 0.2, -0.3],
            ],
            dtype=float,
        )

        def fail_diag(*_args, **_kwargs):
            raise AssertionError("rate extraction should use the Q diagonal view")

        monkeypatch.setattr(scm_module.np, "diag", fail_diag)

        rates, transition_cdfs, deterministic_targets = (
            svc._prepare_branch_history_context(Q, Q.shape[0])
        )

        np.testing.assert_allclose(rates, -Q.diagonal())
        assert len(transition_cdfs) == Q.shape[0]
        assert deterministic_targets == [-1, -1, -1]

    def test_run_single_simulation_uses_precomputed_metadata(
        self, default_args, monkeypatch
    ):
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        svc = StochasticCharacterMap(default_args)
        tip_states = {"A": "x", "B": "x", "C": "y", "D": "y"}
        states = ["x", "y"]
        Q = np.array([[-0.8, 0.8], [0.8, -0.8]])
        pi = np.ones(2) / 2
        parent_map = svc._build_parent_map(tree)
        cond_liks, _ = svc._felsenstein_pruning(tree, tip_states, Q, pi, states)
        metadata = svc._build_simulation_metadata(
            tree, tip_states, states, parent_map
        )

        expected = svc._run_single_simulation(
            tree,
            tip_states,
            Q,
            pi,
            states,
            np.random.default_rng(7),
            cond_liks=cond_liks,
            parent_map=parent_map,
            transition_cache={},
            simulation_metadata=metadata,
        )

        def fail_find_clades(*args, **kwargs):
            raise AssertionError("precomputed metadata should avoid traversal")

        monkeypatch.setattr(tree, "find_clades", fail_find_clades)

        actual = svc._run_single_simulation(
            tree,
            tip_states,
            Q,
            pi,
            states,
            np.random.default_rng(7),
            cond_liks=cond_liks,
            parent_map=parent_map,
            transition_cache={},
            simulation_metadata=metadata,
        )

        assert actual == expected

    def test_run_single_simulation_reuses_branch_history_context(
        self, default_args, monkeypatch
    ):
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        svc = StochasticCharacterMap(default_args)
        tip_states = {"A": "x", "B": "x", "C": "y", "D": "y"}
        states = ["x", "y"]
        Q = np.array([[-0.8, 0.8], [0.8, -0.8]])
        pi = np.ones(2) / 2
        parent_map = svc._build_parent_map(tree)
        cond_liks, _ = svc._felsenstein_pruning(tree, tip_states, Q, pi, states)
        metadata = svc._build_simulation_metadata(
            tree, tip_states, states, parent_map
        )
        sampling_nodes = svc._prepare_ancestral_sampling_nodes(
            metadata[0],
            Q,
            cond_liks,
            len(states),
            {},
        )

        expected = svc._run_single_simulation(
            tree,
            tip_states,
            Q,
            pi,
            states,
            np.random.default_rng(7),
            cond_liks=cond_liks,
            parent_map=parent_map,
            transition_cache={},
            simulation_metadata=metadata,
            sampling_nodes=sampling_nodes,
        )

        branch_history_context = svc._prepare_branch_history_context(Q, len(states))

        def fail_prepare_branch_history_context(*_args, **_kwargs):
            raise AssertionError("branch-history context should be reused")

        monkeypatch.setattr(
            svc,
            "_prepare_branch_history_context",
            fail_prepare_branch_history_context,
        )

        actual = svc._run_single_simulation(
            tree,
            tip_states,
            Q,
            pi,
            states,
            np.random.default_rng(7),
            cond_liks=cond_liks,
            parent_map=parent_map,
            transition_cache={},
            simulation_metadata=metadata,
            sampling_nodes=sampling_nodes,
            branch_history_context=branch_history_context,
        )

        assert actual == expected

    def test_sample_ancestral_states_reuses_transition_cache(
        self, default_args, monkeypatch
    ):
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        svc = StochasticCharacterMap(default_args)
        tip_states = {"A": "x", "B": "x", "C": "y", "D": "y"}
        states = ["x", "y"]
        cond_liks = {
            id(clade): np.ones(2)
            for clade in tree.find_clades(order="preorder")
        }
        parent_map = svc._build_parent_map(tree)
        calls = []

        def fake_matrix_exp(Q, t):
            calls.append(t)
            return np.eye(2)

        monkeypatch.setattr(svc, "_matrix_exp", fake_matrix_exp)
        transition_cache = {}

        node_states = svc._sample_ancestral_states(
            tree,
            tip_states,
            np.array([[-1.0, 1.0], [1.0, -1.0]]),
            np.array([1.0, 0.0]),
            cond_liks,
            states,
            np.random.default_rng(42),
            parent_map=parent_map,
            transition_cache=transition_cache,
        )

        assert calls == [1.0]
        assert set(transition_cache) == {1.0}
        assert node_states[id(tree.root.clades[0])] == 0
        assert node_states[id(tree.root.clades[1])] == 0

    def test_precomputed_ancestral_sampling_preserves_rng_stream(
        self, default_args
    ):
        tree = Phylo.read(StringIO("((A:1,B:2):3,(C:4,D:5):6);"), "newick")
        svc = StochasticCharacterMap(default_args)
        tip_states = {"A": "x", "B": "y", "C": "x", "D": "y"}
        states = ["x", "y", "z"]
        Q = np.array(
            [
                [-0.7, 0.4, 0.3],
                [0.2, -0.5, 0.3],
                [0.1, 0.2, -0.3],
            ],
            dtype=float,
        )
        pi = np.ones(3) / 3
        cond_liks = {
            id(clade): np.ones(3)
            for clade in tree.find_clades(order="preorder")
        }
        for tip in tree.get_terminals():
            idx = 0 if tip_states[tip.name] == "x" else 1
            lik = np.zeros(3)
            lik[idx] = 1.0
            cond_liks[id(tip)] = lik
        parent_map = svc._build_parent_map(tree)
        simulation_nodes, _ = svc._build_simulation_metadata(
            tree,
            tip_states,
            states,
            parent_map,
        )
        transition_cache = {}
        sampling_nodes = svc._prepare_ancestral_sampling_nodes(
            simulation_nodes,
            Q,
            cond_liks,
            len(states),
            transition_cache,
        )

        baseline_rng = np.random.default_rng(2026)
        optimized_rng = np.random.default_rng(2026)
        baseline = svc._sample_ancestral_states(
            tree,
            tip_states,
            Q,
            pi,
            cond_liks,
            states,
            baseline_rng,
            parent_map=parent_map,
            transition_cache=dict(transition_cache),
            simulation_nodes=simulation_nodes,
        )
        optimized = svc._sample_ancestral_states(
            tree,
            tip_states,
            Q,
            pi,
            cond_liks,
            states,
            optimized_rng,
            parent_map=parent_map,
            transition_cache={},
            simulation_nodes=simulation_nodes,
            sampling_nodes=sampling_nodes,
        )

        assert optimized == baseline
        assert optimized_rng.random() == pytest.approx(baseline_rng.random())

    def test_precomputed_ancestral_sampling_avoids_matrix_exponentials(
        self, default_args, monkeypatch
    ):
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        svc = StochasticCharacterMap(default_args)
        tip_states = {"A": "x", "B": "x", "C": "y", "D": "y"}
        states = ["x", "y"]
        Q = np.array([[-1.0, 1.0], [1.0, -1.0]])
        pi = np.ones(2) / 2
        cond_liks = {
            id(clade): np.ones(2)
            for clade in tree.find_clades(order="preorder")
        }
        parent_map = svc._build_parent_map(tree)
        simulation_nodes, _ = svc._build_simulation_metadata(
            tree,
            tip_states,
            states,
            parent_map,
        )
        sampling_nodes = svc._prepare_ancestral_sampling_nodes(
            simulation_nodes,
            Q,
            cond_liks,
            len(states),
            {},
        )

        def fail_matrix_exp(*_args, **_kwargs):
            raise AssertionError("precomputed sampling should not call matrix_exp")

        monkeypatch.setattr(svc, "_matrix_exp", fail_matrix_exp)

        node_states = svc._sample_ancestral_states(
            tree,
            tip_states,
            Q,
            pi,
            cond_liks,
            states,
            np.random.default_rng(42),
            parent_map=parent_map,
            simulation_nodes=simulation_nodes,
            sampling_nodes=sampling_nodes,
        )

        assert set(node_states) == {id(clade) for clade in tree.find_clades()}

    @pytest.mark.parametrize("k", [2, 3, 4])
    def test_prepare_ancestral_sampling_nodes_small_state_paths_match_formula(
        self, default_args, k
    ):
        svc = StochasticCharacterMap(default_args)
        rng = np.random.default_rng(456)
        transition = rng.random((k, k))
        cond_liks = {
            1: rng.random(k),
            2: np.zeros(k),
            3: rng.random(k),
        }
        simulation_nodes = [
            (1, 0, 1.0, None),
            (2, 0, 1.0, None),
            (3, 0, 1.0, 0),
        ]

        expected = []
        uniform_cdf = np.cumsum(np.ones(k) / k)
        for clade_id, parent_id, _branch_length, tip_state_idx in simulation_nodes:
            if tip_state_idx is not None:
                expected.append((clade_id, parent_id, tip_state_idx, None))
                continue

            cdfs_by_parent = np.empty((k, k), dtype=float)
            for parent_state in range(k):
                probs = transition[parent_state, :] * cond_liks[clade_id]
                total = probs.sum()
                if total > 0:
                    cdfs_by_parent[parent_state, :] = np.cumsum(probs / total)
                else:
                    cdfs_by_parent[parent_state, :] = uniform_cdf
            expected.append((clade_id, parent_id, tip_state_idx, cdfs_by_parent))

        observed = svc._prepare_ancestral_sampling_nodes(
            simulation_nodes,
            np.eye(k),
            cond_liks,
            k,
            {1.0: transition},
        )

        for observed_node, expected_node in zip(observed, expected):
            assert observed_node[:3] == expected_node[:3]
            if expected_node[3] is None:
                assert observed_node[3] is None
            else:
                np.testing.assert_allclose(observed_node[3], expected_node[3])

    def test_prepare_ancestral_sampling_nodes_vector_path_matches_row_formula(
        self, default_args
    ):
        svc = StochasticCharacterMap(default_args)
        k = 8
        rng = np.random.default_rng(123)
        transition = rng.random((k, k))
        cond_liks = {
            1: rng.random(k),
            2: np.zeros(k),
            3: rng.random(k),
        }
        simulation_nodes = [
            (1, 0, 1.0, None),
            (2, 0, 1.0, None),
            (3, 0, 1.0, 4),
        ]
        expected = []
        uniform_cdf = np.cumsum(np.ones(k) / k)
        for clade_id, parent_id, _branch_length, tip_state_idx in simulation_nodes:
            if tip_state_idx is not None:
                expected.append((clade_id, parent_id, tip_state_idx, None))
                continue

            cdfs_by_parent = np.empty((k, k), dtype=float)
            for parent_state in range(k):
                probs = transition[parent_state, :] * cond_liks[clade_id]
                total = probs.sum()
                if total > 0:
                    cdfs_by_parent[parent_state, :] = np.cumsum(probs / total)
                else:
                    cdfs_by_parent[parent_state, :] = uniform_cdf
            expected.append((clade_id, parent_id, tip_state_idx, cdfs_by_parent))

        observed = svc._prepare_ancestral_sampling_nodes(
            simulation_nodes,
            np.eye(k),
            cond_liks,
            k,
            {1.0: transition},
        )

        for observed_node, expected_node in zip(observed, expected):
            assert observed_node[:3] == expected_node[:3]
            if expected_node[3] is None:
                assert observed_node[3] is None
            else:
                np.testing.assert_allclose(observed_node[3], expected_node[3])

    def test_prepare_ancestral_sampling_nodes_positive_vector_path_matches_formula(
        self, default_args
    ):
        svc = StochasticCharacterMap(default_args)
        k = 6
        rng = np.random.default_rng(321)
        transition = rng.random((k, k))
        cond_liks = {1: rng.random(k) + 0.1}
        simulation_nodes = [(1, 0, 1.0, None)]

        weighted = transition * cond_liks[1]
        expected_cdfs = np.cumsum(
            weighted / weighted.sum(axis=1)[:, None],
            axis=1,
        )

        observed = svc._prepare_ancestral_sampling_nodes(
            simulation_nodes,
            np.eye(k),
            cond_liks,
            k,
            {1.0: transition},
        )

        assert observed[0][:3] == (1, 0, None)
        np.testing.assert_allclose(observed[0][3], expected_cdfs)

    def test_mapping_produces_valid_histories(self, default_args):
        svc = StochasticCharacterMap(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        tip_states = svc._parse_discrete_trait_file(
            DISCRETE_TRAITS_FILE, "diet", tree_tips
        )
        states = sorted(set(tip_states.values()))
        state_idx = {s: i for i, s in enumerate(states)}
        k = len(states)
        pi = np.ones(k) / k

        Q, _ = svc._fit_q_matrix(tree, tip_states, states, "ER")
        rng = np.random.default_rng(42)
        mapping = svc._run_single_simulation(
            tree, tip_states, Q, pi, states, rng
        )

        assert "node_states" in mapping
        assert "branch_histories" in mapping

        # Check that tip states match
        for tip in tree.get_terminals():
            if tip.name in tip_states:
                expected = state_idx[tip_states[tip.name]]
                assert mapping["node_states"][id(tip)] == expected

        # Check that histories have valid states
        for clade_id, history in mapping["branch_histories"].items():
            assert len(history) >= 1
            for time_val, state in history:
                assert 0 <= state < k
                assert time_val >= 0

    def test_branch_times_sum_correctly(self, default_args):
        svc = StochasticCharacterMap(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        tip_states = svc._parse_discrete_trait_file(
            DISCRETE_TRAITS_FILE, "diet", tree_tips
        )
        states = sorted(set(tip_states.values()))
        k = len(states)
        pi = np.ones(k) / k

        Q, _ = svc._fit_q_matrix(tree, tip_states, states, "ER")
        rng = np.random.default_rng(42)
        mapping = svc._run_single_simulation(
            tree, tip_states, Q, pi, states, rng
        )

        # Sum dwelling times across all branches
        total_time = 0.0
        for clade in tree.find_clades(order="preorder"):
            if clade == tree.root:
                continue
            clade_id = id(clade)
            if clade_id not in mapping["branch_histories"]:
                continue
            t = clade.branch_length if clade.branch_length else 0.0
            history = mapping["branch_histories"][clade_id]
            branch_total = 0.0
            for h_idx in range(len(history)):
                start_t = history[h_idx][0]
                if h_idx + 1 < len(history):
                    end_t = history[h_idx + 1][0]
                else:
                    end_t = t
                branch_total += end_t - start_t
            assert branch_total == pytest.approx(t, abs=1e-6)
            total_time += branch_total

        # Total tree length = 277.2772 (from R)
        assert total_time == pytest.approx(277.2772, abs=0.01)


class TestSummary:
    def test_summarize_simulations_accumulates_with_single_traversal(
        self, default_args, monkeypatch
    ):
        tree = Phylo.read(StringIO("((A:1,B:1):2,C:1);"), "newick")
        internal = tree.root.clades[0]
        internal_id = id(internal)
        branch_ids = [
            id(clade)
            for clade in tree.find_clades(order="preorder")
            if clade != tree.root
        ]
        all_ids = [
            id(clade) for clade in tree.find_clades(order="preorder")
        ]

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("generic traversal should not be used")

        monkeypatch.setattr(tree, "find_clades", fail_find_clades)

        mappings = [
            {
                "branch_histories": {
                    cid: [(0.0, 0)] for cid in branch_ids
                },
                "node_states": {cid: 0 for cid in all_ids},
            },
            {
                "branch_histories": {
                    cid: [(0.0, 1)] for cid in branch_ids
                },
                "node_states": {cid: 1 for cid in all_ids},
            },
        ]
        mappings[0]["branch_histories"][internal_id] = [(0.0, 0), (1.0, 1)]

        svc = StochasticCharacterMap(default_args)
        summary = svc._summarize_simulations(mappings, ["x", "y"], tree)

        assert summary["mean_dwelling_times"] == pytest.approx([2.0, 3.0])
        assert np.allclose(
            summary["mean_transitions"],
            np.array([[0.0, 0.5], [0.0, 0.0]]),
        )
        assert summary["node_posteriors"][id(tree.root)] == pytest.approx([0.5, 0.5])

    def test_summarize_simulations_matches_generic_history_reference(
        self, default_args
    ):
        tree = Phylo.read(StringIO("((A:1,B:1):2,C:1);"), "newick")
        branch_lengths = {
            id(clade): (clade.branch_length if clade.branch_length else 1e-8)
            for clade in tree.find_clades(order="preorder")
            if clade != tree.root
        }
        branch_ids = list(branch_lengths)
        mappings = [
            {
                "branch_histories": {
                    branch_ids[0]: [(0.0, 0)],
                    branch_ids[1]: [(0.0, 1), (0.25, 0)],
                    branch_ids[2]: [(0.0, 1)],
                    branch_ids[3]: [(0.0, 0), (0.5, 1)],
                },
                "node_states": {id(tree.root): 0},
            },
            {
                "branch_histories": {
                    branch_ids[0]: [(0.0, 1)],
                    branch_ids[1]: [(0.0, 0)],
                    branch_ids[2]: [(0.0, 0), (0.4, 1)],
                    branch_ids[3]: [(0.0, 1)],
                },
                "node_states": {id(tree.root): 1},
            },
        ]

        expected_dwelling = np.zeros(2)
        expected_transitions = np.zeros((2, 2))
        for mapping in mappings:
            for clade_id, history in mapping["branch_histories"].items():
                branch_length = branch_lengths[clade_id]
                for h_idx in range(len(history)):
                    state = history[h_idx][1]
                    start_t = history[h_idx][0]
                    if h_idx + 1 < len(history):
                        end_t = history[h_idx + 1][0]
                    else:
                        end_t = branch_length
                    expected_dwelling[state] += end_t - start_t
                for h_idx in range(1, len(history)):
                    from_state = history[h_idx - 1][1]
                    to_state = history[h_idx][1]
                    if from_state != to_state:
                        expected_transitions[from_state, to_state] += 1

        svc = StochasticCharacterMap(default_args)
        summary = svc._summarize_simulations(mappings, ["x", "y"], tree)

        np.testing.assert_allclose(
            summary["mean_dwelling_times"], expected_dwelling / len(mappings)
        )
        np.testing.assert_allclose(
            summary["mean_transitions"], expected_transitions / len(mappings)
        )

    def test_dwelling_times_sum_to_tree_length(self, default_args):
        default_args.nsim = 50
        svc = StochasticCharacterMap(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        tip_states = svc._parse_discrete_trait_file(
            DISCRETE_TRAITS_FILE, "diet", tree_tips
        )
        states = sorted(set(tip_states.values()))
        k = len(states)
        pi = np.ones(k) / k

        Q, _ = svc._fit_q_matrix(tree, tip_states, states, "ER")
        rng = np.random.default_rng(42)
        mappings = []
        for _ in range(50):
            m = svc._run_single_simulation(
                tree, tip_states, Q, pi, states, rng
            )
            mappings.append(m)

        summary = svc._summarize_simulations(mappings, states, tree)
        mean_dwelling = summary["mean_dwelling_times"]

        # Total should equal tree length
        assert mean_dwelling.sum() == pytest.approx(277.2772, abs=0.5)


class TestRun:
    def _stub_expensive_run_tail(self, svc, monkeypatch, captured_tree):
        def fit_q_matrix(tree, *_args, **_kwargs):
            captured_tree["tree"] = tree
            return np.array([[-1.0, 1.0], [1.0, -1.0]]), -1.25

        monkeypatch.setattr(svc, "_fit_q_matrix", fit_q_matrix)
        monkeypatch.setattr(
            svc,
            "_felsenstein_pruning",
            lambda tree, *_args, **_kwargs: ({id(tree.root): np.ones(2)}, None),
        )
        monkeypatch.setattr(svc, "_build_parent_map", lambda *_args, **_kwargs: {})
        monkeypatch.setattr(
            svc,
            "_build_simulation_metadata",
            lambda *_args, **_kwargs: ([],),
        )
        monkeypatch.setattr(
            svc,
            "_prepare_ancestral_sampling_nodes",
            lambda *_args, **_kwargs: [],
        )
        monkeypatch.setattr(
            svc,
            "_prepare_branch_history_context",
            lambda *_args, **_kwargs: {},
        )
        monkeypatch.setattr(
            svc,
            "_summarize_simulations",
            lambda *_args, **_kwargs: {
                "mean_dwelling_times": np.zeros(2),
                "mean_transitions": np.zeros((2, 2)),
            },
        )
        monkeypatch.setattr(svc, "_print_text_output", lambda *_args, **_kwargs: None)

    def test_all_tips_present_uses_read_only_tree_without_copy_or_prune(
        self, default_args, monkeypatch
    ):
        default_args.nsim = 0
        tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )
        tip_states = {"A": "x", "B": "x", "C": "y", "D": "y"}
        svc = StochasticCharacterMap(default_args)
        captured_tree = {}

        monkeypatch.setattr(svc, "read_tree_file_unmodified", lambda: tree)
        monkeypatch.setattr(
            svc,
            "read_tree_file",
            lambda: (_ for _ in ()).throw(
                AssertionError("mutable tree read should not be used")
            ),
        )
        monkeypatch.setattr(
            svc,
            "_parse_discrete_trait_file",
            lambda *_args, **_kwargs: tip_states,
        )
        monkeypatch.setattr(
            svc,
            "_fast_copy",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("tree should not be copied")
            ),
        )
        monkeypatch.setattr(
            svc,
            "_prune_tree_to_tip_states",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("no-prune path should skip pruning")
            ),
        )
        self._stub_expensive_run_tail(svc, monkeypatch, captured_tree)

        svc.run()

        assert captured_tree["tree"] is tree
        assert {tip.name for tip in tree.get_terminals()} == {"A", "B", "C", "D"}

    def test_missing_tip_states_copy_before_pruning(
        self, default_args, monkeypatch
    ):
        default_args.nsim = 0
        tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )
        tip_states = {"A": "x", "B": "x", "C": "y"}
        svc = StochasticCharacterMap(default_args)
        original_fast_copy = svc._fast_copy
        copied_trees = []
        captured_tree = {}

        def copy_spy(tree_to_copy):
            copied_tree = original_fast_copy(tree_to_copy)
            copied_trees.append(copied_tree)
            return copied_tree

        monkeypatch.setattr(svc, "read_tree_file_unmodified", lambda: tree)
        monkeypatch.setattr(
            svc,
            "_parse_discrete_trait_file",
            lambda *_args, **_kwargs: tip_states,
        )
        monkeypatch.setattr(svc, "_fast_copy", copy_spy)
        self._stub_expensive_run_tail(svc, monkeypatch, captured_tree)

        svc.run()

        assert len(copied_trees) == 1
        assert captured_tree["tree"] is copied_trees[0]
        assert {tip.name for tip in tree.get_terminals()} == {"A", "B", "C", "D"}
        assert {tip.name for tip in copied_trees[0].get_terminals()} == {
            "A",
            "B",
            "C",
        }

    def test_text_output(self, default_args, capsys):
        svc = StochasticCharacterMap(default_args)
        svc.run()
        captured = capsys.readouterr()
        assert "Stochastic Character Mapping (SIMMAP)" in captured.out
        assert "Model: ER" in captured.out
        assert "Log-likelihood:" in captured.out
        assert "Mean dwelling times:" in captured.out
        assert "Mean transitions:" in captured.out

    def test_print_text_output_batches_sections(self, default_args):
        svc = StochasticCharacterMap(default_args)
        states = ["herbivore", "carnivore"]
        Q = np.array([[-0.1, 0.1], [0.2, -0.2]])
        summary = {
            "mean_dwelling_times": np.array([3.0, 1.0]),
            "mean_transitions": np.array([[0.0, 2.0], [1.0, 0.0]]),
        }

        with patch("builtins.print") as mocked_print:
            svc._print_text_output(Q, -4.25, states, summary)

        mocked_print.assert_called_once()
        output = mocked_print.call_args.args[0]
        assert "Stochastic Character Mapping (SIMMAP)" in output
        assert "Fitted Q matrix:" in output
        assert "herbivore" in output
        assert "Mean dwelling times:" in output
        assert "herbivore -> carnivore:  2.00" in output
        assert "Total:" in output
        lines = output.splitlines()
        assert lines[7:9] == [
            "herbivore          -0.1000      0.1000",
            "carnivore           0.2000     -0.2000",
        ]
        assert lines[13:15] == [
            "  herbivore          3.00 (75.0%)",
            "  carnivore          1.00 (25.0%)",
        ]
        assert lines[17:19] == [
            "  herbivore -> carnivore:  2.00",
            "  carnivore -> herbivore:  1.00",
        ]

    def test_json_output(self, json_args, capsys):
        svc = StochasticCharacterMap(json_args)
        svc.run()
        captured = capsys.readouterr()
        payload = json.loads(captured.out)
        assert payload["model"] == "ER"
        assert payload["nsim"] == 10
        assert "q_matrix" in payload
        assert "log_likelihood" in payload
        assert "states" in payload
        assert "mean_dwelling_times" in payload
        assert "mean_transitions" in payload
        assert list(payload["q_matrix"]) == payload["states"]
        for row in payload["q_matrix"].values():
            assert list(row) == payload["states"]
            assert all(isinstance(value, float) for value in row.values())
        assert list(payload["mean_dwelling_times"]) == payload["states"]
        assert list(payload["mean_dwelling_proportions"]) == payload["states"]
        expected_transitions = [
            f"{source} -> {target}"
            for source in payload["states"]
            for target in payload["states"]
            if source != target
        ]
        assert list(payload["mean_transitions"]) == expected_transitions


class TestPlot:
    def test_plot_file_created(self, default_args):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(
            suffix=".png", delete=False
        ) as f:
            tmppath = f.name

        try:
            default_args.plot = tmppath
            svc = StochasticCharacterMap(default_args)
            svc.run()
            assert os.path.exists(tmppath)
            assert os.path.getsize(tmppath) > 0
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    @pytest.mark.parametrize("circular", [False, True])
    def test_plot_stochastic_map_uses_direct_tree_traversal(
        self, monkeypatch, default_args, circular
    ):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            default_args.circular = circular
            svc = StochasticCharacterMap(default_args)
            tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
            parent_map = svc._build_parent_map(tree)
            mapping = {"branch_histories": {}, "node_states": {}}

            def fail_traversal(*_args, **_kwargs):
                raise AssertionError("plot setup should use direct traversal")

            monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)
            monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)

            svc._plot_stochastic_map(
                tree, mapping, ["0", "1"], tmppath, parent_map=parent_map
            )

            assert os.path.exists(tmppath)
            assert os.path.getsize(tmppath) > 0
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    def test_iter_preorder_preserves_order_without_reversed(self, default_args):
        class NoReversedList(list):
            def __reversed__(self):
                raise AssertionError("_iter_preorder should push children directly")

        svc = StochasticCharacterMap(default_args)
        root = Clade(name="root")
        left = Clade(name="left")
        middle = Clade(name="middle")
        right = Clade(name="right")
        left.clades = NoReversedList([Clade(name="left_a"), Clade(name="left_b")])
        middle.clades = NoReversedList([Clade(name="middle_a")])
        right.clades = NoReversedList()
        root.clades = NoReversedList([left, middle, right])

        order = [clade.name for clade in svc._iter_preorder(root)]

        assert order == [
            "root",
            "left",
            "left_a",
            "left_b",
            "middle",
            "middle_a",
            "right",
        ]

    @pytest.mark.parametrize("circular", [False, True])
    def test_plot_stochastic_map_skips_redundant_tight_layout(
        self, monkeypatch, default_args, tmp_path, circular
    ):
        pytest.importorskip("matplotlib")
        from matplotlib.figure import Figure

        default_args.circular = circular
        svc = StochasticCharacterMap(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        parent_map = svc._build_parent_map(tree)
        mapping = {"branch_histories": {}, "node_states": {}}

        def fail_tight_layout(*args, **kwargs):
            raise AssertionError("savefig(bbox_inches='tight') handles plot bounds")

        monkeypatch.setattr(Figure, "tight_layout", fail_tight_layout)

        output = str(tmp_path / f"scm_no_tight_layout_{circular}.png")
        svc._plot_stochastic_map(
            tree, mapping, ["0", "1"], output, parent_map=parent_map
        )

        assert os.path.exists(output)
        assert os.path.getsize(output) > 0

    def test_plot_stochastic_map_reuses_preorder_for_node_positions(
        self, monkeypatch, default_args, tmp_path
    ):
        pytest.importorskip("matplotlib")
        import phykit.helpers.plot_config as plot_config

        default_args.circular = False
        svc = StochasticCharacterMap(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        parent_map = svc._build_parent_map(tree)
        expected_ids = [id(clade) for clade in svc._iter_preorder(tree.root)]
        mapping = {"branch_histories": {}, "node_states": {}}
        original_compute_node_positions = plot_config.compute_node_positions
        calls = []

        def assert_preorder_reused(
            tree_arg, parent_map_arg, cladogram=False, preorder_clades=None
        ):
            assert preorder_clades is not None
            calls.append([id(clade) for clade in preorder_clades])
            return original_compute_node_positions(
                tree_arg,
                parent_map_arg,
                cladogram=cladogram,
                preorder_clades=preorder_clades,
            )

        monkeypatch.setattr(
            plot_config, "compute_node_positions", assert_preorder_reused
        )

        output = str(tmp_path / "scm_preorder_positions.png")
        svc._plot_stochastic_map(
            tree, mapping, ["0", "1"], output, parent_map=parent_map
        )

        assert calls == [expected_ids]
        assert os.path.exists(output)
        assert os.path.getsize(output) > 0

    def test_plot_stochastic_map_reuses_clade_lists_for_circular_coords(
        self, monkeypatch, default_args, tmp_path
    ):
        pytest.importorskip("matplotlib")
        import phykit.helpers.circular_layout as circular_layout

        default_args.circular = True
        svc = StochasticCharacterMap(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        parent_map = svc._build_parent_map(tree)
        expected_preorder_ids = [
            id(clade) for clade in svc._iter_preorder(tree.root)
        ]
        expected_terminal_ids = [
            id(clade)
            for clade in svc._iter_preorder(tree.root)
            if not clade.clades
        ]
        mapping = {"branch_histories": {}, "node_states": {}}
        original_compute_circular_coords = circular_layout.compute_circular_coords
        calls = []

        def assert_clade_lists_reused(
            tree_arg,
            node_x_arg,
            parent_map_arg,
            preorder_clades=None,
            terminal_clades=None,
        ):
            assert preorder_clades is not None
            assert terminal_clades is not None
            calls.append(
                (
                    [id(clade) for clade in preorder_clades],
                    [id(clade) for clade in terminal_clades],
                )
            )
            return original_compute_circular_coords(
                tree_arg,
                node_x_arg,
                parent_map_arg,
                preorder_clades=preorder_clades,
                terminal_clades=terminal_clades,
            )

        monkeypatch.setattr(
            circular_layout, "compute_circular_coords", assert_clade_lists_reused
        )

        output = str(tmp_path / "scm_circular_precomputed_coords.png")
        svc._plot_stochastic_map(
            tree, mapping, ["0", "1"], output, parent_map=parent_map
        )

        assert calls == [(expected_preorder_ids, expected_terminal_ids)]
        assert os.path.exists(output)
        assert os.path.getsize(output) > 0

    @pytest.mark.parametrize("circular", [False, True])
    def test_plot_stochastic_map_batches_branch_segments(
        self, monkeypatch, default_args, tmp_path, circular
    ):
        pytest.importorskip("matplotlib")
        import matplotlib.axes
        from matplotlib.collections import LineCollection

        default_args.circular = circular
        svc = StochasticCharacterMap(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        parent_map = svc._build_parent_map(tree)
        clades = list(svc._iter_preorder(tree.root))
        mapping = {
            "branch_histories": {
                id(clade): [(0.0, idx % 2), (0.5, (idx + 1) % 2)]
                for idx, clade in enumerate(clades)
                if clade is not tree.root
            },
            "node_states": {
                id(clade): idx % 2
                for idx, clade in enumerate(clades)
                if clade.clades
            },
        }
        output = str(tmp_path / f"scm_batched_{circular}.png")
        original_add_collection = matplotlib.axes.Axes.add_collection
        line_collections = []

        def fail_plot(*args, **kwargs):
            raise AssertionError("stochastic-map branches should be batched")

        def count_collection(self, collection, *args, **kwargs):
            if isinstance(collection, LineCollection):
                line_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "plot", fail_plot)
        monkeypatch.setattr(
            matplotlib.axes.Axes, "add_collection", count_collection
        )

        svc._plot_stochastic_map(
            tree, mapping, ["0", "1"], output, parent_map=parent_map
        )

        assert len(line_collections) >= 2
        assert os.path.exists(output)
        assert os.path.getsize(output) > 0
