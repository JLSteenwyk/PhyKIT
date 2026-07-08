from argparse import Namespace
from io import StringIO
from pathlib import Path
import pickle as stdlib_pickle
import subprocess
import sys

import pytest
from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin
from Bio.Phylo.Newick import Clade, Tree as NewickTree

from phykit.services.tree import base as base_module
from phykit.services.tree.base import Tree
from phykit.errors import PhykitUserError


def test_tree_base_import_defers_biophylo_and_numpy():
    code = (
        "import sys; "
        "import phykit.services.tree.base; "
        "assert 'typing' not in sys.modules; "
        "assert 'hashlib' not in sys.modules; "
        "assert 'Bio' not in sys.modules; "
        "assert 'Bio.Phylo' not in sys.modules; "
        "assert 'pickle' not in sys.modules; "
        "assert 'numpy' not in sys.modules"
    )

    subprocess.run([sys.executable, "-c", code], check=True)


def test_lazy_pickle_caches_resolved_copy_helpers(monkeypatch):
    lazy_pickle = base_module._LazyPickle()

    def cached_dumps(value, **_kwargs):
        return f"cached:{value}".encode("ascii")

    def cached_loads(value):
        return value.decode("ascii").removeprefix("cached:")

    def uncached_dumps(*_args, **_kwargs):
        raise AssertionError("cached dumps should be reused")

    def uncached_loads(*_args, **_kwargs):
        raise AssertionError("cached loads should be reused")

    monkeypatch.setattr(stdlib_pickle, "dumps", cached_dumps)
    monkeypatch.setattr(stdlib_pickle, "loads", cached_loads)

    protocol = lazy_pickle.HIGHEST_PROTOCOL
    assert lazy_pickle.loads(lazy_pickle.dumps("tree", protocol=protocol)) == "tree"

    monkeypatch.setattr(stdlib_pickle, "dumps", uncached_dumps)
    monkeypatch.setattr(stdlib_pickle, "loads", uncached_loads)

    assert lazy_pickle.loads(lazy_pickle.dumps("tree2", protocol=protocol)) == "tree2"
    assert lazy_pickle.__dict__["dumps"] is cached_dumps
    assert lazy_pickle.__dict__["loads"] is cached_loads
    assert lazy_pickle.__dict__["HIGHEST_PROTOCOL"] == stdlib_pickle.HIGHEST_PROTOCOL


class _Tip:
    def __init__(self, name):
        self.name = name


class _TreeObj:
    def __init__(self):
        self.pruned = []

    def prune(self, taxon):
        self.pruned.append(taxon)

    def get_terminals(self):
        return [_Tip("a"), _Tip("b")]


class _Internal:
    def __init__(self, branch_length):
        self.branch_length = branch_length


class _CalcTree:
    def __init__(self, internals, total_len):
        self._internals = internals
        self._total_len = total_len

    def get_nonterminals(self):
        return self._internals

    def total_branch_length(self):
        return self._total_len


class _TotalBranchFallbackTree:
    def total_branch_length(self):
        return 7.5

    def count_terminals(self):
        return 3


class TestTreeBase:
    def test_init_sets_fields(self):
        service = Tree(
            tree_file_path="tree.tre",
            tree1_file_path="tree1.tre",
            alignment_file_path="aln.fa",
            output_file_path="out.tre",
            reference="ref.tre",
            tip_1="a",
            tip_2="b",
            exclude_gaps=True,
        )
        assert service.tree_file_path == "tree.tre"
        assert service.tree1_file_path == "tree1.tre"
        assert service.alignment_file_path == "aln.fa"
        assert service.output_file_path == "out.tre"
        assert service.reference == "ref.tre"
        assert service.tip_1 == "a"
        assert service.tip_2 == "b"
        assert service.exclude_gaps is True

    def test_get_file_hash_handles_oserror(self, mocker):
        mocker.patch("phykit.services.tree.base.os.stat", side_effect=OSError)
        assert Tree._get_file_hash("/missing") == ""

    def test_get_file_hash_success(self, tmp_path):
        file_path = tmp_path / "t.tre"
        file_path.write_text("(A,B);")
        cache_key = Tree._get_file_hash(str(file_path))
        stat = file_path.stat()
        assert cache_key == f"{file_path}_{stat.st_size}_{stat.st_mtime_ns}"

    def test_calculate_total_branch_length_fast_matches_biophylo(self):
        tree = NewickTree(
            root=Clade(
                branch_length=None,
                clades=[
                    Clade(branch_length=1.5, name="a"),
                    Clade(
                        branch_length=0.0,
                        clades=[
                            Clade(branch_length=2.0, name="b"),
                            Clade(branch_length=-0.5, name="c"),
                        ],
                    ),
                ],
            )
        )

        assert Tree.calculate_total_branch_length_fast(tree) == tree.total_branch_length()

    def test_calculate_total_branch_length_fast_uses_standard_tree_traversal(
        self, monkeypatch
    ):
        tree = NewickTree(
            root=Clade(
                branch_length=None,
                clades=[
                    Clade(branch_length=1.5, name="a"),
                    Clade(branch_length=2.5, name="b"),
                ],
            )
        )

        def fail_total_branch_length(*_args, **_kwargs):
            raise AssertionError("standard trees should use direct traversal")

        monkeypatch.setattr(NewickTree, "total_branch_length", fail_total_branch_length)

        assert Tree.calculate_total_branch_length_fast(tree) == 4.0

    def test_calculate_total_branch_length_fast_skips_empty_child_extension(self):
        class NoIterEmptyList(list):
            def __iter__(self):
                raise AssertionError("terminal child list should not be extended")

        class MinimalClade:
            def __init__(self, branch_length=None, clades=None):
                self.branch_length = branch_length
                self.clades = clades if clades is not None else NoIterEmptyList()

        tree = type(
            "Tree",
            (),
            {
                "root": MinimalClade(
                    branch_length=1.0,
                    clades=[
                        MinimalClade(branch_length=2.0),
                        MinimalClade(branch_length=3.0),
                    ],
                )
            },
        )()

        assert Tree.calculate_total_branch_length_fast(tree) == 6.0

    def test_calculate_total_branch_length_and_terminal_count_fast(self, monkeypatch):
        tree = NewickTree(
            root=Clade(
                branch_length=1.0,
                clades=[
                    Clade(branch_length=2.0, name="a"),
                    Clade(branch_length=3.0, name="b"),
                ],
            )
        )

        def fail_generic_summary(*_args, **_kwargs):
            raise AssertionError("standard trees should use direct traversal")

        monkeypatch.setattr(NewickTree, "total_branch_length", fail_generic_summary)
        monkeypatch.setattr(NewickTree, "count_terminals", fail_generic_summary)

        assert Tree.calculate_total_branch_length_and_terminal_count_fast(tree) == (6.0, 2)

    def test_calculate_terminal_count_fast_matches_biophylo(self):
        tree = NewickTree(
            root=Clade(
                clades=[
                    Clade(name="a"),
                    Clade(
                        clades=[
                            Clade(name="b"),
                            Clade(name="c"),
                        ],
                    ),
                ],
            )
        )

        assert Tree.calculate_terminal_count_fast(tree) == tree.count_terminals()
        assert Tree.calculate_terminal_count_fast(tree.root.clades[1]) == 2

    def test_calculate_terminal_names_fast_preserves_tree_order(self):
        tree = NewickTree(
            root=Clade(
                clades=[
                    Clade(clades=[Clade(name="a"), Clade(name="b")]),
                    Clade(name="c"),
                    Clade(
                        clades=[
                            Clade(name="d"),
                            Clade(name="e"),
                            Clade(name="f"),
                        ],
                    ),
                ],
            )
        )

        assert Tree.calculate_terminal_names_fast(tree) == [
            "a",
            "b",
            "c",
            "d",
            "e",
            "f",
        ]

    def test_calculate_terminal_count_fast_fallback_signal(self):
        tree = _TotalBranchFallbackTree()

        assert Tree.calculate_terminal_count_fast(tree) is None

    def test_total_branch_helpers_fallback_for_tree_like_objects(self):
        tree = _TotalBranchFallbackTree()

        assert Tree.calculate_total_branch_length_fast(tree) == 7.5
        assert Tree.calculate_total_branch_length_and_terminal_count_fast(tree) == (7.5, 3)

    def test_calculate_internal_and_total_branch_length_fast_matches_biophylo(self):
        tree = NewickTree(
            root=Clade(
                branch_length=1.0,
                clades=[
                    Clade(branch_length=2.0, name="a"),
                    Clade(
                        branch_length=3.0,
                        clades=[
                            Clade(branch_length=4.0, name="b"),
                            Clade(branch_length=None, name="c"),
                        ],
                    ),
                ],
            )
        )
        expected_internal = sum(
            clade.branch_length
            for clade in tree.get_nonterminals()
            if clade.branch_length is not None
        )

        assert Tree.calculate_internal_and_total_branch_length_fast(tree) == (
            expected_internal,
            tree.total_branch_length(),
        )

    def test_calculate_internal_and_total_branch_length_fast_fallback(self):
        tree = _CalcTree([_Internal(2.0), _Internal(None), _Internal(1.0)], 6.0)

        assert Tree.calculate_internal_and_total_branch_length_fast(tree) == (3.0, 6.0)

    def test_calculate_terminal_root_distances_fast(self):
        tree = NewickTree(
            root=Clade(
                branch_length=5.0,
                clades=[
                    Clade(branch_length=1.0, name="a"),
                    Clade(
                        branch_length=2.0,
                        clades=[
                            Clade(branch_length=3.0, name="b"),
                            Clade(branch_length=None, name="c"),
                        ],
                    ),
                ],
            )
        )

        distances = Tree.calculate_terminal_root_distances_fast(tree)

        assert distances.tolist() == [1.0, 5.0, 2.0]

    def test_calculate_terminal_root_distances_fast_preserves_mixed_child_order(self):
        tree = NewickTree(
            root=Clade(
                branch_length=None,
                clades=[
                    Clade(branch_length=1.0, name="a"),
                    Clade(
                        branch_length=2.0,
                        clades=[
                            Clade(branch_length=3.0, name="b"),
                            Clade(branch_length=4.0, name="c"),
                            Clade(branch_length=5.0, name="d"),
                        ],
                    ),
                    Clade(branch_length=6.0, name="e"),
                ],
            )
        )

        distances = Tree.calculate_terminal_root_distances_fast(tree)

        assert distances.tolist() == [1.0, 5.0, 6.0, 7.0, 6.0]

    def test_validate_tree_uses_direct_traversal_for_branch_lengths(self, monkeypatch):
        service = Tree()
        tree = Phylo.read(StringIO("((a:1,b:1):1,c:1);"), "newick")

        def fail_traversal(*_args, **_kwargs):
            raise AssertionError("generic tree traversal should not be used")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)
        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)

        service.validate_tree(tree, min_tips=3, require_branch_lengths=True)

    def test_validate_tree_assigns_default_lengths_with_direct_traversal(self, monkeypatch):
        service = Tree()
        tree = NewickTree(
            root=Clade(
                branch_length=None,
                clades=[
                    Clade(branch_length=None, name="a"),
                    Clade(branch_length=1.0, name="b"),
                    Clade(branch_length=None, name="c"),
                ],
            )
        )

        def fail_traversal(*_args, **_kwargs):
            raise AssertionError("generic tree traversal should not be used")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)
        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)

        service.validate_tree(
            tree,
            min_tips=3,
            assign_default_branch_length=1e-8,
        )

        assert tree.root.branch_length is None
        assert [clade.branch_length for clade in tree.root.clades] == [
            1e-8,
            1.0,
            1e-8,
        ]

    def test_validate_tree_handles_mixed_child_counts_with_direct_traversal(
        self, monkeypatch
    ):
        service = Tree()
        tree = NewickTree(
            root=Clade(
                branch_length=None,
                clades=[
                    Clade(branch_length=None, name="a"),
                    Clade(
                        branch_length=None,
                        clades=[
                            Clade(branch_length=1.0, name="b"),
                            Clade(branch_length=None, name="c"),
                            Clade(branch_length=2.0, name="d"),
                        ],
                    ),
                    Clade(branch_length=3.0, name="e"),
                ],
            )
        )

        def fail_traversal(*_args, **_kwargs):
            raise AssertionError("generic tree traversal should not be used")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)
        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)

        service.validate_tree(
            tree,
            min_tips=5,
            assign_default_branch_length=1e-8,
        )

        assert tree.root.branch_length is None
        assert tree.root.clades[0].branch_length == 1e-8
        assert tree.root.clades[1].branch_length == 1e-8
        assert [clade.branch_length for clade in tree.root.clades[1].clades] == [
            1.0,
            1e-8,
            2.0,
        ]

    def test_validate_tree_standard_min_tips_stops_after_threshold(self):
        service = Tree()

        class UnvisitedClade:
            @property
            def clades(self):
                raise AssertionError("validation should stop after min_tips")

        tree = NewickTree(
            root=Clade(
                clades=[
                    UnvisitedClade(),
                    Clade(name="a"),
                    Clade(name="b"),
                    Clade(name="c"),
                ],
            )
        )

        service.validate_tree(tree, min_tips=3)

    def test_validate_tree_fallback_stops_after_minimum_tip_count(self):
        service = Tree()

        class NonstandardTree:
            @property
            def root(self):
                raise AttributeError

            def get_terminals(self):
                for idx in range(10):
                    if idx > 2:
                        raise AssertionError("fallback should stop at min_tips")
                    yield object()

        service.validate_tree(NonstandardTree(), min_tips=3)

    def test_validate_tree_fallback_rejects_too_few_tips(self):
        service = Tree()

        class NonstandardTree:
            @property
            def root(self):
                raise AttributeError

            def get_terminals(self):
                yield object()
                yield object()

        with pytest.raises(PhykitUserError):
            service.validate_tree(NonstandardTree(), min_tips=3)

    def test_read_tree_file_success_returns_deepcopy(self, mocker):
        service = Tree(tree_file_path="x.tre")
        source_tree = {"name": "tree"}
        mocker.patch.object(Tree, "_get_file_hash", return_value="h")
        mocker.patch.object(Tree, "_cached_tree_read", return_value=source_tree)
        out = service.read_tree_file()
        assert out == source_tree
        assert out is not source_tree

    def test_read_tree_file_unmodified_returns_cached_tree(self, mocker):
        service = Tree(tree_file_path="x.tre")
        source_tree = {"name": "tree"}
        mocker.patch.object(Tree, "_get_file_hash", return_value="h")
        mocker.patch.object(Tree, "_cached_tree_read", return_value=source_tree)

        out = service.read_tree_file_unmodified()

        assert out is source_tree

    def test_read_tree1_file_not_found_exits(self, capsys):
        service = Tree(tree1_file_path="/missing/tree1.tre")
        with pytest.raises(PhykitUserError) as exc:
            service.read_tree1_file()
        assert exc.value.code == 2
        assert "corresponds to no such file or directory" in exc.value.messages[0]

    def test_read_tree1_file_success_returns_deepcopy(self, mocker):
        service = Tree(tree1_file_path="x1.tre")
        source_tree = {"name": "tree1"}
        mocker.patch.object(Tree, "_get_file_hash", return_value="h")
        mocker.patch.object(Tree, "_cached_tree_read", return_value=source_tree)
        out = service.read_tree1_file()
        assert out == source_tree
        assert out is not source_tree

    def test_read_reference_tree_file_not_found_exits(self, capsys):
        service = Tree(reference="/missing/ref.tre")
        with pytest.raises(PhykitUserError) as exc:
            service.read_reference_tree_file()
        assert exc.value.code == 2
        assert "corresponds to no such file or directory" in exc.value.messages[0]

    def test_read_reference_tree_file_success_returns_deepcopy(self, mocker):
        service = Tree(reference="ref.tre")
        source_tree = {"name": "ref"}
        mocker.patch.object(Tree, "_get_file_hash", return_value="h")
        mocker.patch.object(Tree, "_cached_tree_read", return_value=source_tree)
        out = service.read_reference_tree_file()
        assert out == source_tree
        assert out is not source_tree

    def test_read_tree_with_error_uses_requested_attr_name(self, capsys):
        service = Tree(tree_file_path="/missing/tree.tre")
        with pytest.raises(PhykitUserError) as exc:
            service._read_tree_with_error("/missing/tree.tre", "tree_file_path")
        assert exc.value.code == 2
        assert "/missing/tree.tre corresponds to no such file or directory." in exc.value.messages[0]

    def test_write_tree_file_delegates_to_phylo(self, mocker):
        service = Tree()
        mocked_write = mocker.patch("Bio.Phylo.write", return_value=1)
        tree = object()
        res = service.write_tree_file(tree, "out.tre")
        assert res == 1
        mocked_write.assert_called_once_with(tree, "out.tre", "newick")

    def test_get_tip_names_from_tree(self):
        service = Tree()
        assert service.get_tip_names_from_tree(_TreeObj()) == ["a", "b"]

    def test_get_tip_names_from_tree_uses_direct_terminal_names(self, monkeypatch):
        service = Tree()
        tree = Phylo.read(StringIO("((a:1,b:1):1,(c:1,d:1):1);"), "newick")

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("standard tree should use direct terminal traversal")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        assert service.get_tip_names_from_tree(tree) == ["a", "b", "c", "d"]

    def test_tips_to_prune_for_ordered_mapping_uses_ordered_tail(self):
        class MappingLike(dict):
            def __contains__(self, key):
                raise AssertionError("ordered path should not use membership")

        values = MappingLike({"a": 1.0, "b": 2.0, "c": 3.0})

        assert Tree._tips_to_prune_for_ordered_mapping(
            ["a", "b", "c", "d", "e"],
            values,
            min_ordered_size=0,
        ) == ["d", "e"]

    def test_tips_to_prune_for_ordered_mapping_iterates_prefix(self):
        class NoIntegerIndexList(list):
            def __getitem__(self, key):
                if isinstance(key, int):
                    raise AssertionError(
                        "ordered prefix comparison should not index tree tips"
                    )
                return super().__getitem__(key)

        values = {"a": 1.0, "b": 2.0, "c": 3.0}

        assert Tree._tips_to_prune_for_ordered_mapping(
            NoIntegerIndexList(["a", "b", "c", "d", "e"]),
            values,
            min_ordered_size=0,
        ) == ["d", "e"]

    def test_tips_to_prune_for_ordered_mapping_falls_back_for_interleaved_tips(self):
        values = {"a": 1.0, "b": 2.0, "c": 3.0}

        assert Tree._tips_to_prune_for_ordered_mapping(
            ["a", "d", "b", "c"],
            values,
            min_ordered_size=0,
        ) == ["d"]

    def test_tips_to_prune_for_ordered_mapping_keeps_small_input_membership_path(self):
        values = {"a": 1.0, "b": 2.0, "c": 3.0}

        assert Tree._tips_to_prune_for_ordered_mapping(
            ["a", "b", "c", "d"],
            values,
        ) == ["d"]

    def test_tips_to_prune_for_ordered_names_uses_ordered_tail(self):
        class OrderedNames(list):
            def __contains__(self, key):
                raise AssertionError("ordered path should not use membership")

        values = OrderedNames(["a", "b", "c"])

        assert Tree._tips_to_prune_for_ordered_names(
            ["a", "b", "c", "d", "e"],
            values,
            min_ordered_size=0,
        ) == ["d", "e"]

    def test_tips_to_prune_for_ordered_names_falls_back_for_interleaved_tips(self):
        assert Tree._tips_to_prune_for_ordered_names(
            ["a", "d", "b", "c"],
            ["a", "b", "c"],
            min_ordered_size=0,
        ) == ["d"]

    def test_tips_to_prune_for_ordered_names_keeps_small_input_set_path(self):
        class OrderedNames(list):
            def __contains__(self, key):
                raise AssertionError("small fallback should use a set")

        values = OrderedNames(["a", "b", "c"])

        assert Tree._tips_to_prune_for_ordered_names(
            ["a", "b", "c", "d"],
            values,
        ) == ["d"]

    def test_get_first_tip_name_from_tree_uses_direct_leftmost_terminal(
        self, monkeypatch
    ):
        service = Tree()
        tree = Phylo.read(StringIO("(((a:1,b:1):1,c:1):1,(d:1,e:1):1);"), "newick")

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("standard tree should use direct first-tip traversal")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        assert service.get_first_tip_name_from_tree(tree) == "a"

    def test_calculate_first_terminal_name_fast_fallback_signal(self):
        class NonstandardClade:
            name = "x"
            clades = ("not", "a", "list")

        assert Tree.calculate_first_terminal_name_fast(NonstandardClade()) is None

    def test_calculate_first_terminal_name_fast_fallback_for_inner_nonstandard_clade(self):
        class InnerNonstandardClade:
            name = "ignored"
            clades = ("not", "a", "list")

        class RootClade:
            clades = [InnerNonstandardClade()]

        assert Tree.calculate_first_terminal_name_fast(RootClade()) is None

    def test_calculate_terminal_names_fast_preserves_tree_and_clade_order(self):
        tree = Phylo.read(StringIO("((a:1,b:1):1,(c:1,d:1):1);"), "newick")

        assert Tree.calculate_terminal_names_fast(tree) == [
            tip.name for tip in tree.get_terminals()
        ]
        assert Tree.calculate_terminal_names_fast(tree.root.clades[1]) == ["c", "d"]

    def test_calculate_terminal_clades_fast_preserves_tree_order(self):
        tree = Phylo.read(StringIO("((a:1,b:1):1,(c:1,d:1):1);"), "newick")

        terminals = Tree.calculate_terminal_clades_fast(tree)

        assert [terminal.name for terminal in terminals] == ["a", "b", "c", "d"]
        assert terminals == tree.get_terminals()

    def test_calculate_terminal_clades_fast_preserves_mixed_child_order(self):
        tree = NewickTree(
            root=Clade(
                clades=[
                    Clade(clades=[Clade(name="a"), Clade(name="b")]),
                    Clade(name="c"),
                    Clade(
                        clades=[
                            Clade(name="d"),
                            Clade(name="e"),
                            Clade(name="f"),
                        ],
                    ),
                ],
            )
        )
        expected = tree.get_terminals()

        terminals = Tree.calculate_terminal_clades_fast(tree)

        assert [terminal.name for terminal in terminals] == [
            "a",
            "b",
            "c",
            "d",
            "e",
            "f",
        ]
        assert terminals == expected

    def test_shared_tips(self):
        service = Tree()
        shared = service.shared_tips(["a", "b"], ["b", "c"])
        assert shared == ["b"]

    def test_shared_tips_preserves_unique_set_semantics(self):
        service = Tree()
        shared = service.shared_tips(["a", "b", "b"], ["b", "b", "c"])
        assert shared == ["b"]

    def test_shared_tips_no_overlap_exits(self, capsys):
        service = Tree()
        with pytest.raises(PhykitUserError) as exc:
            service.shared_tips(["a"], ["b"])
        assert exc.value.code == 2
        assert exc.value.messages == ["no common tips"]

    def test_prune_tree_using_taxa_list(self):
        service = Tree()
        tree = _TreeObj()
        res = service.prune_tree_using_taxa_list(tree, ["a", "b"])
        assert res is tree
        assert tree.pruned == ["a", "b"]

    def test_prune_tree_using_taxa_list_empty_list_is_noop(self):
        service = Tree()
        tree = _TreeObj()
        res = service.prune_tree_using_taxa_list(tree, [])
        assert res is tree
        assert tree.pruned == []

    def test_prune_tree_using_taxa_list_uses_terminal_objects(self, monkeypatch):
        service = Tree()
        tree = Phylo.read(StringIO("((a:1,b:1):1,(c:1,d:1):1);"), "newick")
        original_find_any = TreeMixin.find_any

        def fail_string_find_any(self, *args, **kwargs):
            if args and args[0] == "a":
                raise AssertionError("string prune lookup should be avoided")
            return original_find_any(self, *args, **kwargs)

        monkeypatch.setattr(TreeMixin, "find_any", fail_string_find_any)

        res = service.prune_tree_using_taxa_list(tree, ["a"])

        assert res is tree
        assert {tip.name for tip in tree.get_terminals()} == {"b", "c", "d"}

    def test_prune_tree_using_taxa_list_uses_direct_terminal_index(self, monkeypatch):
        service = Tree()
        tree = Phylo.read(StringIO("((a:1,b:1):1,(c:1,d:1):1);"), "newick")

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("generic terminal traversal should not be used")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        res = service.prune_tree_using_taxa_list(tree, ["a"])

        assert res is tree
        assert Tree.calculate_terminal_names_fast(tree) == ["b", "c", "d"]

    def test_terminal_targets_and_count_fast_handles_mixed_child_counts(self):
        tree = NewickTree(
            root=Clade(
                clades=[
                    Clade(name="a"),
                    Clade(
                        clades=[
                            Clade(name="b"),
                            Clade(name="c"),
                            Clade(name="d"),
                        ],
                    ),
                    Clade(name="e"),
                ],
            )
        )

        targets, terminal_count = Tree._terminal_targets_and_count_fast(
            tree,
            ["d", "a", "e"],
        )

        assert terminal_count == 5
        assert [target.name for target in targets] == ["d", "a", "e"]

    def test_terminal_by_name_fast_handles_mixed_child_counts(self):
        tree = NewickTree(
            root=Clade(
                clades=[
                    Clade(name="a"),
                    Clade(
                        clades=[
                            Clade(name="b"),
                            Clade(name="c"),
                            Clade(name="d"),
                        ],
                    ),
                    Clade(name="e"),
                ],
            )
        )

        terminal_by_name = Tree._terminal_by_name_fast(tree)

        assert set(terminal_by_name) == {"a", "b", "c", "d", "e"}
        assert terminal_by_name["c"].name == "c"

    def test_prune_tree_using_taxa_list_only_indexes_requested_terminals(
        self, monkeypatch
    ):
        service = Tree()
        tree = Phylo.read(StringIO("((a:1,b:1):1,(c:1,d:1):1);"), "newick")
        tree.root.clades[1].clades[0].name = ["unhashable", "non-target"]

        def fail_generic_path(*_args, **_kwargs):
            raise AssertionError("selected target lookup should use batch pruning")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_generic_path)
        monkeypatch.setattr(TreeMixin, "prune", fail_generic_path)

        res = service.prune_tree_using_taxa_list(tree, ["a"])

        assert res is tree
        names = Tree.calculate_terminal_names_fast(tree)
        assert "a" not in names
        assert "b" in names
        assert ["unhashable", "non-target"] in names

    def test_prune_tree_using_taxa_list_batches_standard_tree_pruning(
        self, monkeypatch
    ):
        service = Tree()
        tree = Phylo.read(StringIO("((a:1,b:2):3,(c:4,d:5):6);"), "newick")

        def fail_prune(*_args, **_kwargs):
            raise AssertionError("per-target prune should not be used")

        monkeypatch.setattr(TreeMixin, "prune", fail_prune)

        res = service.prune_tree_using_taxa_list(tree, ["a", "c"])

        assert res is tree
        assert Tree.calculate_terminal_names_fast(tree) == ["b", "d"]
        assert [tip.branch_length for tip in tree.get_terminals()] == [5.0, 11.0]

    def test_prune_tree_using_taxa_list_batch_matches_biophylo_prune(self):
        service = Tree()
        newick = "((a:1,b:2,c:3):4,(d:5,e:6):7,f:8);"
        expected = Phylo.read(StringIO(newick), "newick")
        actual = Phylo.read(StringIO(newick), "newick")
        taxa_to_prune = ["a", "b", "d", "f"]
        expected_by_name = {tip.name: tip for tip in expected.get_terminals()}
        for taxon in taxa_to_prune:
            expected.prune(expected_by_name[taxon])

        service.prune_tree_using_taxa_list(actual, taxa_to_prune)

        expected_out = StringIO()
        actual_out = StringIO()
        Phylo.write(expected, expected_out, "newick")
        Phylo.write(actual, actual_out, "newick")
        assert actual_out.getvalue() == expected_out.getvalue()

    def test_batch_prune_helper_matches_biophylo_without_generic_paths(
        self, monkeypatch
    ):
        newick = "((a:1,b:2,c:3):4,(d:5,e:6):7,f:8);"
        expected = Phylo.read(StringIO(newick), "newick")
        actual = Phylo.read(StringIO(newick), "newick")
        taxa_to_prune = {"a", "b", "d", "f"}
        expected_by_name = {tip.name: tip for tip in expected.get_terminals()}
        for taxon in taxa_to_prune:
            expected.prune(expected_by_name[taxon])
        target_ids = {
            id(tip) for tip in actual.get_terminals() if tip.name in taxa_to_prune
        }

        def fail_generic_path(*_args, **_kwargs):
            raise AssertionError("batch helper should not use generic tree methods")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_generic_path)
        monkeypatch.setattr(TreeMixin, "prune", fail_generic_path)

        assert Tree._prune_terminal_objects_batch_standard_tree(actual, target_ids)

        expected_out = StringIO()
        actual_out = StringIO()
        Phylo.write(expected, expected_out, "newick")
        Phylo.write(actual, actual_out, "newick")
        assert actual_out.getvalue() == expected_out.getvalue()

    def test_prune_tree_using_taxa_list_preserves_all_tip_prune_error(self):
        service = Tree()
        tree = Phylo.read(StringIO("(a:1,b:1);"), "newick")

        with pytest.raises(ValueError):
            service.prune_tree_using_taxa_list(tree, ["a", "b"])

    def test_calculate_pairwise_tip_distances_fast_uses_direct_traversal(self, monkeypatch):
        service = Tree()
        tree = Phylo.read(StringIO("(a:1,(b:0.5,c:0.5):0.5);"), "newick")

        def fail_get_path(*args, **kwargs):
            raise AssertionError("get_path should not be called")

        def fail_generic_traversal(*args, **kwargs):
            raise AssertionError("generic tree traversal should not be used")

        monkeypatch.setattr(TreeMixin, "get_path", fail_get_path)
        monkeypatch.setattr(TreeMixin, "depths", fail_generic_traversal)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_generic_traversal)

        combos, distances = service.calculate_pairwise_tip_distances_fast(
            tree, ["a", "b", "c"]
        )

        assert combos == [("a", "b"), ("a", "c"), ("b", "c")]
        assert distances == pytest.approx([2.0, 2.0, 1.0])

        combos, distances = service.calculate_pairwise_tip_distances_fast(
            tree,
            ["a", "b", "c"],
            include_combos=False,
        )

        assert combos is None
        assert distances == pytest.approx([2.0, 2.0, 1.0])

    def test_calculate_pairwise_tip_distances_fast_handles_mixed_child_counts(
        self, monkeypatch
    ):
        service = Tree()
        tree = NewickTree(
            root=Clade(
                clades=[
                    Clade(branch_length=1.0, name="a"),
                    Clade(
                        branch_length=2.0,
                        clades=[
                            Clade(branch_length=3.0, name="b"),
                            Clade(branch_length=4.0, name="c"),
                            Clade(branch_length=5.0, name="d"),
                        ],
                    ),
                    Clade(branch_length=6.0, name="e"),
                ],
            )
        )

        def fail_generic_traversal(*_args, **_kwargs):
            raise AssertionError("generic tree traversal should not be used")

        monkeypatch.setattr(TreeMixin, "get_path", fail_generic_traversal)
        monkeypatch.setattr(TreeMixin, "depths", fail_generic_traversal)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_generic_traversal)

        combos, distances = service.calculate_pairwise_tip_distances_fast(
            tree,
            ["e", "b", "d", "a"],
        )

        assert combos == [
            ("e", "b"),
            ("e", "d"),
            ("e", "a"),
            ("b", "d"),
            ("b", "a"),
            ("d", "a"),
        ]
        assert distances == pytest.approx([11.0, 13.0, 7.0, 8.0, 6.0, 8.0])

    def test_calculate_pairwise_tip_distances_fast_deep_tree_uses_lca_index(self, monkeypatch):
        service = Tree()
        root = Clade()
        current = root
        tips = []
        for idx in range(7):
            tip_name = f"tip{idx}"
            tips.append(tip_name)
            next_node = Clade(branch_length=1.0)
            current.clades = [
                Clade(branch_length=1.0, name=tip_name),
                next_node,
            ]
            current = next_node
        tips.append("tip7")
        current.clades = [Clade(branch_length=1.0, name="tip7")]
        tree = NewickTree(root=root)

        def fail_path_helper(*_args, **_kwargs):
            raise AssertionError("deep all-pairs distances should use LCA index")

        monkeypatch.setattr(Tree, "_PAIRWISE_LCA_DEPTH_THRESHOLD", 4)
        monkeypatch.setattr(
            Tree,
            "_pairwise_tip_distances_from_paths",
            staticmethod(fail_path_helper),
        )

        combos, distances = service.calculate_pairwise_tip_distances_fast(tree, tips)
        expected_combos = [
            (tips[i], tips[j])
            for i in range(len(tips) - 1)
            for j in range(i + 1, len(tips))
        ]
        expected_distances = [
            tree.distance(tip_a, tip_b)
            for tip_a, tip_b in expected_combos
        ]

        assert combos == expected_combos
        assert distances == pytest.approx(expected_distances)

        combos, distances = service.calculate_pairwise_tip_distances_fast(
            tree,
            tips,
            include_combos=False,
        )

        assert combos is None
        assert distances == pytest.approx(expected_distances)

    def test_calculate_treeness_print_value(self, capsys):
        service = Tree()
        tree = _CalcTree([_Internal(2.0), _Internal(None), _Internal(1.0)], total_len=6.0)
        val = service.calculate_treeness(tree=tree, print_value=True)
        assert val == 0.5
        out, _ = capsys.readouterr()
        assert out.strip() == "0.5"

    def test_calculate_treeness_reads_tree_when_not_supplied(self, mocker):
        service = Tree(tree_file_path="x.tre")
        tree = _CalcTree([_Internal(1.0)], total_len=2.0)
        mocker.patch.object(Tree, "read_tree_file", return_value=tree)
        assert service.calculate_treeness() == 0.5

    def test_calculate_treeness_zero_length(self, capsys):
        service = Tree()
        tree = _CalcTree([_Internal(1.0)], total_len=0.0)
        val = service.calculate_treeness(tree=tree)
        assert val is None
        out, _ = capsys.readouterr()
        assert "Invalid tree. Tree should contain branch lengths" in out

    def test_calculate_treeness_print_broken_pipe(self, mocker):
        service = Tree()
        tree = _CalcTree([_Internal(1.0)], total_len=2.0)
        mocker.patch("builtins.print", side_effect=BrokenPipeError)
        val = service.calculate_treeness(tree=tree, print_value=True)
        assert val is None

    def test_calculate_treeness_zero_length_print_broken_pipe(self, mocker):
        service = Tree()
        tree = _CalcTree([_Internal(1.0)], total_len=0.0)
        mocker.patch("builtins.print", side_effect=BrokenPipeError)
        val = service.calculate_treeness(tree=tree)
        assert val is None

    def test_get_gap_chars(self):
        assert Tree.get_gap_chars(is_protein=True) == ["-", "?", "*", "X", "x"]
        assert Tree.get_gap_chars(is_protein=False) == ["-", "?", "*", "X", "x", "N", "n"]
