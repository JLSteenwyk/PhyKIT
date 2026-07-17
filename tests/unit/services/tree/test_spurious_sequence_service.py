from argparse import Namespace
from io import StringIO
import subprocess
import sys

import pytest
from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin

import phykit.services.tree.spurious_sequence as spurious_sequence_module
from phykit.services.tree.spurious_sequence import SpuriousSequence


def test_module_import_does_not_import_statistics_or_json():
    code = """
import sys
import phykit.services.tree.spurious_sequence as module
assert callable(module.print_json)
assert hasattr(module.np, "__getattr__")
assert "typing" not in sys.modules
assert "statistics" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "numpy" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_lazy_numpy_caches_resolved_attributes():
    lazy_np = spurious_sequence_module._LazyNumpy()

    first_asarray = lazy_np.asarray
    second_asarray = lazy_np.asarray

    assert first_asarray is second_asarray
    assert lazy_np.__dict__["asarray"] is first_asarray
    assert lazy_np._module is not None


@pytest.fixture
def args():
    return Namespace(tree="/some/path/to/file.tre", factor=None)


class _Terminal:
    def __init__(self, name, branch_length):
        self.name = name
        self.branch_length = branch_length


class _Tree:
    def __init__(self, terminals):
        self._terminals = terminals

    def get_terminals(self):
        return self._terminals


class _NoSpuriousBranchMap:
    def values(self):
        return [1.0, 2.0]

    def items(self):
        raise AssertionError("no row scan needed when no branch meets threshold")


class TestSpuriousSequence:
    def test_init_sets_expected_attrs(self, args):
        service = SpuriousSequence(args)
        assert service.tree_file_path == args.tree
        assert service.factor == 20
        assert service.json_output is False
        assert service.method == "median-factor"
        assert service.alpha == 0.05
        assert service.max_remove is None
        assert service.tree_list is False
        assert service.per_species is False

    @pytest.mark.parametrize(
        "values, message",
        [
            (
                {"factor": 2, "method": "diameter-impact"},
                "--factor cannot be used",
            ),
            (
                {"method": "diameter-impact", "alpha": 0},
                "--alpha must be",
            ),
            (
                {"method": "diameter-impact", "max_remove": 0},
                "--max-remove must be",
            ),
            (
                {"method": "diameter-impact", "per_species": True},
                "--per-species requires",
            ),
            (
                {"method": "median-factor", "tree_list": True},
                "require --method diameter-impact",
            ),
            (
                {"method": "median-factor", "alpha": 0.1},
                "require --method diameter-impact",
            ),
        ],
    )
    def test_process_args_rejects_incompatible_diameter_options(
        self,
        values,
        message,
    ):
        arguments = {
            "tree": "tree.tre",
            "factor": None,
            "method": "median-factor",
            "alpha": 0.05,
            "max_remove": None,
            "tree_list": False,
            "per_species": False,
            "json": False,
        }
        arguments.update(values)

        with pytest.raises(SystemExit) as error:
            SpuriousSequence(Namespace(**arguments))

        assert message in error.value.messages[0]

    def test_get_branch_lengths_and_names_uses_terminal_branches_only(self, args):
        service = SpuriousSequence(args)
        tree = _Tree(
            [
                _Terminal("a", 1.0),
                _Terminal("b", None),
                _Terminal("c", 3.0),
            ]
        )
        branch_lengths, name_map = service.get_branch_lengths_and_their_names(tree)
        assert branch_lengths == [1.0, 3.0]
        assert name_map == {"a": 1.0, "c": 3.0}

    def test_get_branch_lengths_fast_path_does_not_call_get_terminals(
        self, args, monkeypatch
    ):
        service = SpuriousSequence(args)
        tree = Phylo.read(StringIO("((a:1,b:2):3,c:4);"), "newick")

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("standard tree should use direct traversal")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        branch_lengths, name_map = service.get_branch_lengths_and_their_names(tree)

        assert branch_lengths == [1.0, 2.0, 4.0]
        assert name_map == {"a": 1.0, "b": 2.0, "c": 4.0}

    def test_get_branch_lengths_standard_tree_collects_in_one_pass(
        self, args, monkeypatch
    ):
        service = SpuriousSequence(args)
        tree = Phylo.read(StringIO("((a:1,b:2):3,c:4);"), "newick")

        def fail_terminal_list(_tree):
            raise AssertionError(
                "standard tree should collect branch lengths directly"
            )

        monkeypatch.setattr(
            SpuriousSequence,
            "_iter_terminal_clades",
            staticmethod(fail_terminal_list),
        )

        branch_lengths, name_map = service.get_branch_lengths_and_their_names(tree)

        assert branch_lengths == [1.0, 2.0, 4.0]
        assert name_map == {"a": 1.0, "b": 2.0, "c": 4.0}

    def test_iter_terminal_clades_preserves_order_with_mixed_child_counts(
        self, monkeypatch
    ):
        tree = Phylo.read(
            StringIO("(a:1,(b:2,c:3):4,(d:5,e:6,f:7):8);"),
            "newick",
        )

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("standard tree should use direct traversal")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        terminals = SpuriousSequence._iter_terminal_clades(tree)

        assert [terminal.name for terminal in terminals] == [
            "a",
            "b",
            "c",
            "d",
            "e",
            "f",
        ]

    def test_identify_spurious_sequence(self, args):
        service = SpuriousSequence(args)
        tree = _Tree([_Terminal("a", 1.0), _Terminal("b", 3.0), _Terminal("c", 5.0)])
        name_map, threshold, median = service.identify_spurious_sequence(tree, factor=2)
        assert median == 3.0
        assert threshold == 6.0
        assert name_map["c"] == 5.0

    def test_median_branch_length_matches_statistics_median(self):
        assert SpuriousSequence._median_branch_length([3.0, 1.0, 5.0]) == 3.0
        assert SpuriousSequence._median_branch_length([4.0, 1.0, 3.0, 2.0]) == 2.5

    def test_median_branch_length_uses_partition_for_large_inputs(self, monkeypatch):
        calls = []

        class _FakeArray:
            def __init__(self, values):
                self.values = sorted(values)

            def partition(self, kth):
                calls.append(kth)

            def __getitem__(self, index):
                return self.values[index]

        class _FakeNumpy:
            def asarray(self, values, dtype=float):
                assert dtype is float
                return _FakeArray(values)

        monkeypatch.setattr(spurious_sequence_module, "_MEDIAN_NUMPY_THRESHOLD", 4)
        monkeypatch.setattr(spurious_sequence_module, "np", _FakeNumpy())

        assert SpuriousSequence._median_branch_length([4.0, 1.0, 3.0, 2.0]) == 2.5
        assert calls == [(1, 2)]

    def test_median_branch_length_default_threshold_uses_partition(self, monkeypatch):
        calls = []

        class _FakeArray:
            def __init__(self, values):
                self.values = sorted(values)

            def partition(self, kth):
                calls.append(kth)

            def __getitem__(self, index):
                return self.values[index]

        class _FakeNumpy:
            def asarray(self, values, dtype=float):
                assert dtype is float
                return _FakeArray(values)

        monkeypatch.setattr(spurious_sequence_module, "np", _FakeNumpy())

        values = list(range(spurious_sequence_module._MEDIAN_NUMPY_THRESHOLD))
        assert SpuriousSequence._median_branch_length(values) == 511.5
        assert calls == [(511, 512)]

    def test_has_spurious_sequence_falls_back_for_nan_max(self):
        assert SpuriousSequence._has_spurious_sequence(
            {"nan": float("nan"), "hit": 100.0},
            10.0,
        )
        assert not SpuriousSequence._has_spurious_sequence(
            {"nan": float("nan"), "miss": 1.0},
            10.0,
        )

    def test_run_prints_none_when_no_spurious_sequences(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", factor=20, json=False)
        service = SpuriousSequence(args)
        mocker.patch.object(
            SpuriousSequence,
            "read_tree_file_unmodified",
            return_value=_Tree([]),
        )
        mocker.patch.object(
            SpuriousSequence,
            "identify_spurious_sequence",
            return_value=({"a": 1.0, "b": 2.0}, 100.0, 1.5),
        )
        mocked_print = mocker.patch("builtins.print")
        service.run()
        mocked_print.assert_called_once_with("None")

    def test_run_skips_text_row_scan_when_no_spurious_sequences(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", factor=20, json=False)
        service = SpuriousSequence(args)
        mocker.patch.object(
            SpuriousSequence,
            "read_tree_file_unmodified",
            return_value=_Tree([]),
        )
        mocker.patch.object(
            SpuriousSequence,
            "identify_spurious_sequence",
            return_value=(_NoSpuriousBranchMap(), 10.0, 1.5),
        )
        mocked_print = mocker.patch("builtins.print")

        service.run()

        mocked_print.assert_called_once_with("None")

    def test_run_prints_spurious_sequences(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", factor=2, json=False)
        service = SpuriousSequence(args)
        mocker.patch.object(
            SpuriousSequence,
            "read_tree_file_unmodified",
            return_value=_Tree([]),
        )
        mocker.patch.object(
            SpuriousSequence,
            "identify_spurious_sequence",
            return_value=({"a": 10.0, "b": 8.0, "c": 1.0}, 5.0, 2.5),
        )
        mocked_print = mocker.patch("builtins.print")
        service.run()
        mocked_print.assert_called_once_with("a\t10.0\t5.0\t2.5\nb\t8.0\t5.0\t2.5")
        assert "c" not in mocked_print.call_args.args[0]

    def test_run_json_payload(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", factor=2, json=True)
        service = SpuriousSequence(args)
        mocker.patch.object(
            SpuriousSequence,
            "read_tree_file_unmodified",
            return_value=_Tree([]),
        )
        mocker.patch.object(
            SpuriousSequence,
            "identify_spurious_sequence",
            return_value=({"a": 10.12345, "b": 1.0}, 5.0, 2.5),
        )
        mocked_json = mocker.patch("phykit.services.tree.spurious_sequence.print_json")
        service.run()
        payload = mocked_json.call_args.args[0]
        assert payload["rows"] == payload["spurious_sequences"]
        assert payload["rows"][0] == {
            "taxon": "a",
            "branch_length": 10.1235,
            "threshold": 5.0,
            "median": 2.5,
        }

    def test_run_skips_json_row_scan_when_no_spurious_sequences(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", factor=20, json=True)
        service = SpuriousSequence(args)
        mocker.patch.object(
            SpuriousSequence,
            "read_tree_file_unmodified",
            return_value=_Tree([]),
        )
        mocker.patch.object(
            SpuriousSequence,
            "identify_spurious_sequence",
            return_value=(_NoSpuriousBranchMap(), 10.0, 1.5),
        )
        mocked_json = mocker.patch("phykit.services.tree.spurious_sequence.print_json")

        service.run()

        assert mocked_json.call_args.args[0] == {
            "rows": [],
            "spurious_sequences": [],
        }

    def test_run_uses_unmodified_tree_read(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", factor=2, json=False)
        tree = _Tree([])
        service = SpuriousSequence(args)
        read_tree = mocker.patch.object(
            service,
            "read_tree_file_unmodified",
            return_value=tree,
        )
        mocker.patch.object(
            service,
            "identify_spurious_sequence",
            return_value=({"a": 1.0}, 10.0, 1.0),
        )
        mocked_print = mocker.patch("builtins.print")

        service.run()

        read_tree.assert_called_once_with()
        mocked_print.assert_called_once_with("None")

    def test_default_max_remove_matches_treeshrink_rule(self):
        assert SpuriousSequence._default_max_remove(4) == 1
        assert SpuriousSequence._default_max_remove(100) == 25
        assert SpuriousSequence._default_max_remove(10_000) == 500

    def test_diameter_shrink_trajectory_selects_high_impact_tip(self, args):
        tree = Phylo.read(
            StringIO(
                "(outlier:100,"
                + ",".join(f"tip_{index}:1" for index in range(1, 20))
                + ");"
            ),
            "newick",
        )
        service = SpuriousSequence(args)
        adjacency, names, _lengths = service._tree_graph(tree)

        trajectory, initial_diameter = service._diameter_shrink_trajectory(
            adjacency,
            names,
            max_remove=5,
        )

        assert initial_diameter == 101.0
        assert trajectory[0]["taxon"] == "outlier"
        assert trajectory[0]["diameter_after"] == 2.0
        assert trajectory[0]["signature"] == pytest.approx(3.9219733363)

    def test_per_gene_calibration_flags_extreme_signature(self, args):
        service = SpuriousSequence(args)
        analysis = {
            "scores": {
                "outlier": 3.9219733363,
                **{f"tip_{index}": 0.0 for index in range(1, 20)},
            },
            "p_values": {},
        }

        service._assign_per_gene_p_values(analysis)

        assert analysis["p_values"]["outlier"] < 0.05
        assert analysis["p_values"]["tip_1"] > 0.05

    def test_collection_kde_is_conservative_for_small_groups(self):
        assert SpuriousSequence._kde_survival_probabilities([0.0, 0.0, 4.0]) == [
            1.0,
            1.0,
            1.0,
        ]

    def test_collection_kde_flags_repeated_species_outlier(self):
        probabilities = SpuriousSequence._kde_survival_probabilities(
            [4.0, 0.0, 0.0, 0.0, 0.0]
        )

        assert probabilities[0] < 0.05
        assert all(probability > 0.05 for probability in probabilities[1:])

    def test_read_tree_list_resolves_relative_paths_and_ignores_comments(
        self,
        tmp_path,
    ):
        list_path = tmp_path / "trees.txt"
        list_path.write_text("# gene trees\n\ngene1.tre\nsub/gene2.tre\n")

        paths = SpuriousSequence._read_tree_list(str(list_path))

        assert paths == [
            str(tmp_path / "gene1.tre"),
            str(tmp_path / "sub" / "gene2.tre"),
        ]
