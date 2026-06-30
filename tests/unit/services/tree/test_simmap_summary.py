from argparse import Namespace
import json
import subprocess
import sys

import pytest
import numpy as np
from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin
from io import StringIO
from unittest.mock import patch

from phykit.services.tree.simmap_summary import SimmapSummary


def test_module_import_does_not_import_numpy():
    code = """
import sys
import phykit.services.tree.simmap_summary as module
assert hasattr(module.np, "__getattr__")
assert callable(module.print_json)
assert "typing" not in sys.modules
assert "json" not in sys.modules
assert "numpy" not in sys.modules
assert "matplotlib" not in sys.modules
assert "matplotlib.pyplot" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


TREE = "tests/sample_files/tree_simple.tre"
TRAITS = "tests/sample_files/tree_simple_discrete_traits.tsv"


def _make_args(**overrides):
    defaults = dict(
        tree=TREE, trait_data=TRAITS, trait="diet",
        model="ER", nsim=20, seed=42, plot=None, csv=None, json=False,
        fig_width=None, fig_height=None, dpi=150, no_title=False,
        title=None, legend_position=None, ylabel_fontsize=None,
        xlabel_fontsize=None, title_fontsize=None, axis_fontsize=None,
        colors=None, ladderize=False, cladogram=False, circular=False,
        color_file=None,
    )
    defaults.update(overrides)
    return Namespace(**defaults)


class TestSimmapSummary:
    def test_branch_label_uses_cached_tip_names(self, monkeypatch):
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        internal = next(
            clade
            for clade in tree.find_clades(order="preorder")
            if not clade.is_terminal() and clade != tree.root
        )
        clade_tip_names = SimmapSummary._collect_clade_tip_names(tree)

        def fail_get_terminals():
            raise AssertionError("cached label path should not call get_terminals")

        monkeypatch.setattr(internal, "get_terminals", fail_get_terminals)

        assert SimmapSummary._branch_label(
            internal,
            tree,
            parent_map=None,
            clade_tip_names=clade_tip_names,
        ) == "(A, B)"

    def test_collect_clade_tip_names_uses_direct_traversal(self, monkeypatch):
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("generic traversal should not be used")

        monkeypatch.setattr(tree, "find_clades", fail_find_clades)

        clade_tip_names = SimmapSummary._collect_clade_tip_names(tree)

        assert clade_tip_names[id(tree.root)] == ("A", "B", "C", "D")
        assert clade_tip_names[id(tree.root.clades[0])] == ("A", "B")
        assert clade_tip_names[id(tree.root.clades[1])] == ("C", "D")

    def test_collect_clade_tip_names_matches_generic_reference(self, monkeypatch):
        tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);"),
            "newick",
        )
        expected = {}
        for clade in tree.find_clades(order="postorder"):
            if clade.is_terminal():
                expected[id(clade)] = (clade.name,)
            else:
                tips = []
                for child in clade.clades:
                    tips.extend(expected[id(child)])
                expected[id(clade)] = tuple(sorted(tips))

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("standard tree cache should use direct traversal")

        monkeypatch.setattr(tree, "find_clades", fail_find_clades)

        observed = SimmapSummary._collect_clade_tip_names(tree)

        assert observed == expected

    def test_iter_preorder_binary_children_use_indexed_path(self):
        from Bio.Phylo.BaseTree import Clade

        class IndexedOnlyList(list):
            def __reversed__(self):
                raise AssertionError("preorder traversal should index children")

        left = Clade(
            name="left",
            clades=IndexedOnlyList([Clade(name="A"), Clade(name="B")]),
        )
        right = Clade(
            name="right",
            clades=IndexedOnlyList([Clade(name="C"), Clade(name="D")]),
        )
        root = Clade(name="root", clades=IndexedOnlyList([left, right]))

        clades = list(SimmapSummary._iter_preorder(root))

        assert clades == [
            root,
            left,
            left.clades[0],
            left.clades[1],
            right,
            right.clades[0],
            right.clades[1],
        ]

    def test_summarize_per_branch_accumulates_with_single_traversal(
        self, monkeypatch
    ):
        svc = SimmapSummary(_make_args())
        tree = Phylo.read(StringIO("((A:1,B:1):2,C:1);"), "newick")
        internal = tree.root.clades[0]
        internal_id = id(internal)
        branch_ids = [
            id(clade)
            for clade in tree.find_clades(order="preorder")
            if clade != tree.root
        ]

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("standard branch summary should use direct traversal")

        monkeypatch.setattr(tree, "find_clades", fail_find_clades)
        mappings = [
            {
                "branch_histories": {
                    cid: [(0.0, 0)] for cid in branch_ids
                },
                "node_states": {},
            },
            {
                "branch_histories": {
                    cid: [(0.0, 1)] for cid in branch_ids
                },
                "node_states": {},
            },
        ]
        mappings[0]["branch_histories"][internal_id] = [(0.0, 0), (1.0, 1)]

        result = svc._summarize_per_branch(mappings, ["x", "y"], tree, {})

        assert result[internal_id]["mean_dwelling"] == pytest.approx([0.5, 1.5])
        assert result[internal_id]["proportions"] == pytest.approx([0.25, 0.75])
        assert result[internal_id]["mean_transitions"] == pytest.approx(0.5)

    def test_summarize_per_branch_matches_generic_history_reference(self):
        svc = SimmapSummary(_make_args())
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
                    branch_ids[1]: [(0.0, 1), (0.2, 0)],
                    branch_ids[2]: [(0.0, 0), (0.4, 1), (0.8, 0)],
                    branch_ids[3]: [(0.0, 1)],
                },
                "node_states": {},
            },
            {
                "branch_histories": {
                    branch_ids[0]: [(0.0, 1)],
                    branch_ids[1]: [(0.0, 0)],
                    branch_ids[2]: [(0.0, 1)],
                    branch_ids[3]: [(0.0, 0), (0.5, 1)],
                },
                "node_states": {},
            },
        ]

        expected_dwelling = {
            cid: np.zeros(2) for cid in branch_lengths
        }
        expected_transitions = {cid: 0.0 for cid in branch_lengths}
        for mapping in mappings:
            for cid, history in mapping["branch_histories"].items():
                branch_length = branch_lengths[cid]
                for h_idx in range(len(history)):
                    state = history[h_idx][1]
                    start_t = history[h_idx][0]
                    if h_idx + 1 < len(history):
                        end_t = history[h_idx + 1][0]
                    else:
                        end_t = branch_length
                    expected_dwelling[cid][state] += end_t - start_t
                for h_idx in range(1, len(history)):
                    if history[h_idx][1] != history[h_idx - 1][1]:
                        expected_transitions[cid] += 1

        result = svc._summarize_per_branch(mappings, ["x", "y"], tree, {})

        for cid in branch_lengths:
            np.testing.assert_allclose(
                result[cid]["mean_dwelling"],
                expected_dwelling[cid] / len(mappings),
            )
            assert result[cid]["mean_transitions"] == pytest.approx(
                expected_transitions[cid] / len(mappings)
            )

    def test_summarize_node_posteriors_uses_direct_traversal(self, monkeypatch):
        svc = SimmapSummary(_make_args())
        tree = Phylo.read(StringIO("((A:1,B:1):2,C:1);"), "newick")
        internal_ids = [
            id(tree.root),
            id(tree.root.clades[0]),
        ]

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("generic traversal should not be used")

        monkeypatch.setattr(tree, "find_clades", fail_find_clades)

        mappings = [
            {"node_states": {internal_ids[0]: 0, internal_ids[1]: 1}},
            {"node_states": {internal_ids[0]: 1, internal_ids[1]: 1}},
        ]

        posteriors = svc._summarize_node_posteriors(
            mappings, ["x", "y"], tree
        )

        np.testing.assert_allclose(posteriors[internal_ids[0]], [0.5, 0.5])
        np.testing.assert_allclose(posteriors[internal_ids[1]], [0.0, 1.0])

    def test_summarize_node_posteriors_handles_mixed_child_counts(self, monkeypatch):
        svc = SimmapSummary(_make_args())
        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1):2,(D:1,E:1,F:1):3);"),
            "newick",
        )
        root = tree.root
        binary = root.clades[1]
        trifurcation = root.clades[2]
        internal_ids = [id(root), id(binary), id(trifurcation)]

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("generic traversal should not be used")

        monkeypatch.setattr(tree, "find_clades", fail_find_clades)

        mappings = [
            {"node_states": {internal_ids[0]: 0, internal_ids[1]: 1}},
            {"node_states": {internal_ids[0]: 1, internal_ids[2]: 0}},
        ]

        posteriors = svc._summarize_node_posteriors(
            mappings,
            ["x", "y"],
            tree,
        )

        assert list(posteriors) == internal_ids
        np.testing.assert_allclose(posteriors[internal_ids[0]], [0.5, 0.5])
        np.testing.assert_allclose(posteriors[internal_ids[1]], [0.0, 0.5])
        np.testing.assert_allclose(posteriors[internal_ids[2]], [0.5, 0.0])

    def test_basic_run(self, capsys):
        """Runs without error and prints summary."""
        svc = SimmapSummary(_make_args())
        svc.run()
        captured = capsys.readouterr()
        assert "SIMMAP Summary" in captured.out
        assert "Q matrix:" in captured.out
        assert "Per-branch dwelling" in captured.out
        assert "Node posterior" in captured.out

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
            lambda *_args, **_kwargs: ([], []),
        )
        monkeypatch.setattr(
            svc,
            "_summarize_per_branch",
            lambda *_args, **_kwargs: {},
        )
        monkeypatch.setattr(
            svc,
            "_summarize_node_posteriors",
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
        monkeypatch.setattr(svc, "_print_text", lambda *_args, **_kwargs: None)

    def test_all_tips_present_uses_read_only_tree_without_copy_or_prune(
        self, monkeypatch
    ):
        tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )
        tip_states = {"A": "x", "B": "x", "C": "y", "D": "y"}
        svc = SimmapSummary(_make_args(nsim=0))
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

    def test_missing_tip_states_copy_before_pruning(self, monkeypatch):
        tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )
        tip_states = {"A": "x", "B": "x", "C": "y"}
        svc = SimmapSummary(_make_args(nsim=0))
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

    @patch("builtins.print")
    def test_print_text_batches_branch_and_node_rows(self, mocked_print, monkeypatch):
        svc = SimmapSummary(_make_args())
        tree = Phylo.read(StringIO("((A:1,B:1):2,C:1);"), "newick")
        states = ["x", "y"]
        Q = np.array([[-0.2, 0.2], [0.1, -0.1]])
        branch_summary = {}
        node_posteriors = {}
        for clade in svc._iter_preorder(tree.root):
            if clade != tree.root:
                branch_summary[id(clade)] = {
                    "branch_length": clade.branch_length or 0.0,
                    "proportions": np.array([0.25, 0.75]),
                    "mean_transitions": 0.5,
                }
            if clade.clades:
                node_posteriors[id(clade)] = np.array([0.4, 0.6])
        tree_summary = {
            "mean_transitions": np.array([[0.0, 1.25], [0.5, 0.0]])
        }

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("text output should reuse direct preorder")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        svc._print_text(
            Q,
            -3.5,
            states,
            branch_summary,
            node_posteriors,
            tree_summary,
            tree,
            {},
        )

        mocked_print.assert_called_once()
        output = mocked_print.call_args.args[0]
        assert "SIMMAP Summary (20 stochastic maps)" in output
        assert "Q matrix:" in output
        assert "x -> y:  1.25" in output
        assert "Per-branch dwelling time proportions:" in output
        assert "Node posterior probabilities:" in output
        assert "(A, B)" in output
        assert "C" in output
        lines = output.splitlines()
        assert lines[7:9] == [
            "x                  -0.2000      0.2000",
            "y                   0.1000     -0.1000",
        ]
        assert lines[11:13] == [
            "  x -> y:  1.25",
            "  y -> x:  0.50",
        ]
        assert lines[18:22] == [
            "(A, B)                                2.0000     0.2500     0.7500       0.50",
            "A                                     1.0000     0.2500     0.7500       0.50",
            "B                                     1.0000     0.2500     0.7500       0.50",
            "C                                     1.0000     0.2500     0.7500       0.50",
        ]
        assert lines[26:28] == [
            "(A, B, C)                               0.4000     0.6000",
            "(A, B)                                  0.4000     0.6000",
        ]

    def test_print_json_output_uses_single_direct_preorder(
        self, monkeypatch
    ):
        svc = SimmapSummary(_make_args(json=True))
        tree = Phylo.read(StringIO("((A:1,B:1):2,C:1);"), "newick")
        states = ["x", "y"]
        Q = np.array([[-0.2, 0.2], [0.1, -0.1]])
        branch_summary = {}
        node_posteriors = {}
        for clade in svc._iter_preorder(tree.root):
            if clade != tree.root:
                branch_summary[id(clade)] = {
                    "branch_length": clade.branch_length or 0.0,
                    "proportions": np.array([0.25, 0.75]),
                    "mean_transitions": 0.5,
                }
            if clade.clades:
                node_posteriors[id(clade)] = np.array([0.4, 0.6])
        tree_summary = {
            "mean_transitions": np.array([[0.0, 1.25], [0.5, 0.0]])
        }
        captured = {}

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("JSON output should reuse direct preorder")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)
        monkeypatch.setattr(
            "phykit.services.tree.simmap_summary.print_json",
            lambda payload, **_kwargs: captured.setdefault("payload", payload),
        )

        svc._print_json_output(
            Q,
            -3.5,
            states,
            branch_summary,
            node_posteriors,
            tree_summary,
            tree,
            {},
        )

        payload = captured["payload"]
        assert [branch["branch"] for branch in payload["branches"]] == [
            "(A, B)",
            "A",
            "B",
            "C",
        ]
        assert [node["node"] for node in payload["node_posteriors"]] == [
            "(A, B, C)",
            "(A, B)",
        ]

    def test_log_likelihood_matches_r(self, capsys):
        """Log-likelihood matches phytools::fitMk (ER) closely."""
        svc = SimmapSummary(_make_args())
        svc.run()
        captured = capsys.readouterr()
        # R gives -8.7889; PhyKIT should be close
        assert "-8.78" in captured.out or "-8.79" in captured.out

    def test_all_states_in_output(self, capsys):
        """All states appear in output."""
        svc = SimmapSummary(_make_args())
        svc.run()
        captured = capsys.readouterr()
        assert "carnivore" in captured.out
        assert "herbivore" in captured.out
        assert "omnivore" in captured.out

    def test_per_branch_proportions_sum_to_one(self, capsys):
        """Dwelling proportions on each branch should sum to ~1."""
        svc = SimmapSummary(_make_args(json=True))
        svc.run()
        captured = capsys.readouterr()
        for line in captured.out.strip().split("\n"):
            if line.startswith("{"):
                payload = json.loads(line)
                break
        for branch in payload["branches"]:
            total = sum(branch["dwelling_proportions"].values())
            assert total == pytest.approx(1.0, abs=0.01), (
                f"Branch {branch['branch']} proportions sum to {total}"
            )

    def test_node_posteriors_sum_to_one(self, capsys):
        """Posteriors at each node should sum to 1."""
        svc = SimmapSummary(_make_args(json=True))
        svc.run()
        captured = capsys.readouterr()
        for line in captured.out.strip().split("\n"):
            if line.startswith("{"):
                payload = json.loads(line)
                break
        for node in payload["node_posteriors"]:
            total = sum(node["posteriors"].values())
            assert total == pytest.approx(1.0, abs=0.01)

    def test_json_output_structure(self, capsys):
        """JSON has correct top-level keys."""
        svc = SimmapSummary(_make_args(json=True))
        svc.run()
        captured = capsys.readouterr()
        for line in captured.out.strip().split("\n"):
            if line.startswith("{"):
                payload = json.loads(line)
                break
        assert payload["model"] == "ER"
        assert payload["nsim"] == 20
        assert "q_matrix" in payload
        assert "branches" in payload
        assert "node_posteriors" in payload
        assert "total_expected_transitions" in payload
        assert len(payload["states"]) == 3
        assert list(payload["q_matrix"]) == payload["states"]
        for row in payload["q_matrix"].values():
            assert list(row) == payload["states"]
            assert all(isinstance(value, float) for value in row.values())

    def test_csv_output(self, tmp_path):
        """CSV file is created with expected content."""
        csv_path = str(tmp_path / "simmap.csv")
        svc = SimmapSummary(_make_args(csv=csv_path))
        svc.run()
        with open(csv_path) as f:
            content = f.read()
        assert "branch,branch_length" in content
        assert "prop_carnivore" in content
        assert "posterior_carnivore" in content

    def test_write_csv_uses_single_direct_preorder(self, monkeypatch, tmp_path):
        svc = SimmapSummary(_make_args())
        tree = Phylo.read(StringIO("((A:1,B:1):2,C:1);"), "newick")
        states = ["x", "y", "z"]
        branch_summary = {}
        node_posteriors = {}
        for clade in svc._iter_preorder(tree.root):
            if clade != tree.root:
                branch_summary[id(clade)] = {
                    "branch_length": clade.branch_length or 0.0,
                    "proportions": np.array([0.25, 0.5, 0.25]),
                    "mean_transitions": 0.5,
                }
            if clade.clades:
                node_posteriors[id(clade)] = np.array([0.4, 0.35, 0.25])

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("CSV output should reuse direct preorder")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        csv_path = tmp_path / "simmap_direct.csv"
        svc._write_csv(
            branch_summary,
            node_posteriors,
            states,
            tree,
            {},
            str(csv_path),
        )

        lines = csv_path.read_text().splitlines()
        assert lines[:5] == [
            "branch,branch_length,prop_x,prop_y,prop_z,expected_transitions",
            "(A; B),2.000000,0.2500,0.5000,0.2500,0.5000",
            "A,1.000000,0.2500,0.5000,0.2500,0.5000",
            "B,1.000000,0.2500,0.5000,0.2500,0.5000",
            "C,1.000000,0.2500,0.5000,0.2500,0.5000",
        ]
        assert lines[-2:] == [
            "(A; B; C),0.4000,0.3500,0.2500",
            "(A; B),0.4000,0.3500,0.2500",
        ]

    def test_plot_output(self, tmp_path):
        """Plot file is created."""
        plot_path = str(tmp_path / "simmap_pie.png")
        svc = SimmapSummary(_make_args(plot=plot_path))
        svc.run()
        assert (tmp_path / "simmap_pie.png").exists()

    def test_plot_posterior_pie_uses_direct_tree_traversal(
        self, monkeypatch, tmp_path
    ):
        svc = SimmapSummary(_make_args())
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        parent_map = svc._build_parent_map(tree)
        clades = list(svc._iter_preorder(tree.root))
        states = ["x", "y"]
        branch_summary = {
            id(clade): {"proportions": np.array([0.7, 0.3])}
            for clade in clades
            if clade != tree.root
        }
        node_posteriors = {
            id(clade): np.array([0.4, 0.6])
            for clade in clades
            if clade.clades and clade != tree.root
        }

        def fail_traversal(*args, **kwargs):
            raise AssertionError("plot setup should use direct traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)

        output_path = tmp_path / "simmap_direct.png"
        svc._plot_posterior_pie(
            tree, node_posteriors, states, branch_summary,
            parent_map, str(output_path),
        )

        assert output_path.exists()

    def test_plot_posterior_pie_reuses_preorder_for_node_positions(
        self, monkeypatch, tmp_path
    ):
        pytest.importorskip("matplotlib")
        import phykit.helpers.plot_config as plot_config

        svc = SimmapSummary(_make_args())
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        parent_map = svc._build_parent_map(tree)
        clades = list(svc._iter_preorder(tree.root))
        states = ["x", "y"]
        branch_summary = {
            id(clade): {"proportions": np.array([0.7, 0.3])}
            for clade in clades
            if clade != tree.root
        }
        node_posteriors = {
            id(clade): np.array([0.4, 0.6])
            for clade in clades
            if clade.clades and clade != tree.root
        }
        expected_preorder_ids = [id(clade) for clade in clades]
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

        output_path = tmp_path / "simmap_preorder_positions.png"
        svc._plot_posterior_pie(
            tree, node_posteriors, states, branch_summary,
            parent_map, str(output_path),
        )

        assert calls == [expected_preorder_ids]
        assert output_path.exists()

    def test_plot_posterior_pie_batches_base_branches(
        self, monkeypatch, tmp_path
    ):
        pytest.importorskip("matplotlib")
        import matplotlib.axes
        from matplotlib.collections import LineCollection

        svc = SimmapSummary(_make_args())
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        parent_map = svc._build_parent_map(tree)
        clades = list(svc._iter_preorder(tree.root))
        states = ["x", "y"]
        branch_summary = {
            id(clade): {"proportions": np.array([0.7, 0.3])}
            for clade in clades
            if clade != tree.root
        }
        node_posteriors = {
            id(clade): np.array([0.4, 0.6])
            for clade in clades
            if clade.clades and clade != tree.root
        }
        output_path = tmp_path / "simmap_batched.png"
        original_add_collection = matplotlib.axes.Axes.add_collection
        line_collections = []

        def fail_plot(*args, **kwargs):
            raise AssertionError("SIMMAP summary branches should be batched")

        def count_collection(self, collection, *args, **kwargs):
            if isinstance(collection, LineCollection):
                line_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "plot", fail_plot)
        monkeypatch.setattr(
            matplotlib.axes.Axes, "add_collection", count_collection
        )

        svc._plot_posterior_pie(
            tree, node_posteriors, states, branch_summary,
            parent_map, str(output_path),
        )

        assert len(line_collections) >= 2
        assert output_path.exists()

    def test_reproducible_with_seed(self, capsys):
        """Same seed gives same results."""
        svc1 = SimmapSummary(_make_args(json=True, seed=123))
        svc1.run()
        out1 = capsys.readouterr().out

        svc2 = SimmapSummary(_make_args(json=True, seed=123))
        svc2.run()
        out2 = capsys.readouterr().out

        j1 = j2 = None
        for line in out1.strip().split("\n"):
            if line.startswith("{"):
                j1 = json.loads(line)
        for line in out2.strip().split("\n"):
            if line.startswith("{"):
                j2 = json.loads(line)
        assert j1["total_expected_transitions"] == j2["total_expected_transitions"]
        assert j1["branches"][0] == j2["branches"][0]

    def test_ard_model(self, capsys):
        """ARD model runs without error."""
        svc = SimmapSummary(_make_args(model="ARD", nsim=10))
        svc.run()
        captured = capsys.readouterr()
        assert "SIMMAP Summary" in captured.out

    def test_expected_transitions_positive(self, capsys):
        """Total expected transitions should be positive."""
        svc = SimmapSummary(_make_args(json=True))
        svc.run()
        captured = capsys.readouterr()
        for line in captured.out.strip().split("\n"):
            if line.startswith("{"):
                payload = json.loads(line)
                break
        assert payload["total_expected_transitions"] > 0

    def test_tip_branches_have_correct_dominant_state(self, capsys):
        """Terminal branches should be dominated by the observed tip state."""
        svc = SimmapSummary(_make_args(json=True, nsim=50))
        svc.run()
        captured = capsys.readouterr()
        for line in captured.out.strip().split("\n"):
            if line.startswith("{"):
                payload = json.loads(line)
                break
        # raccoon is carnivore - its branch should be mostly carnivore
        for b in payload["branches"]:
            if b["branch"] == "raccoon":
                assert b["dwelling_proportions"]["carnivore"] > 0.4
                break
