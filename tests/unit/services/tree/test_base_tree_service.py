from argparse import Namespace
from pathlib import Path

import pytest

from phykit.services.tree.base import Tree


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
        digest = Tree._get_file_hash(str(file_path))
        assert isinstance(digest, str)
        assert len(digest) == 32

    def test_read_tree_file_success_returns_deepcopy(self, mocker):
        service = Tree(tree_file_path="x.tre")
        source_tree = {"name": "tree"}
        mocker.patch.object(Tree, "_get_file_hash", return_value="h")
        mocker.patch.object(Tree, "_cached_tree_read", return_value=source_tree)
        out = service.read_tree_file()
        assert out == source_tree
        assert out is not source_tree

    def test_read_tree1_file_not_found_exits(self, capsys):
        service = Tree(tree1_file_path="/missing/tree1.tre")
        with pytest.raises(SystemExit) as exc:
            service.read_tree1_file()
        assert exc.value.code == 2
        out, _ = capsys.readouterr()
        assert "corresponds to no such file or directory" in out

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
        with pytest.raises(SystemExit) as exc:
            service.read_reference_tree_file()
        assert exc.value.code == 2
        out, _ = capsys.readouterr()
        assert "corresponds to no such file or directory" in out

    def test_read_reference_tree_file_success_returns_deepcopy(self, mocker):
        service = Tree(reference="ref.tre")
        source_tree = {"name": "ref"}
        mocker.patch.object(Tree, "_get_file_hash", return_value="h")
        mocker.patch.object(Tree, "_cached_tree_read", return_value=source_tree)
        out = service.read_reference_tree_file()
        assert out == source_tree
        assert out is not source_tree

    def test_write_tree_file_delegates_to_phylo(self, mocker):
        service = Tree()
        mocked_write = mocker.patch("phykit.services.tree.base.Phylo.write", return_value=1)
        tree = object()
        res = service.write_tree_file(tree, "out.tre")
        assert res == 1
        mocked_write.assert_called_once_with(tree, "out.tre", "newick")

    def test_get_tip_names_from_tree(self):
        service = Tree()
        assert service.get_tip_names_from_tree(_TreeObj()) == ["a", "b"]

    def test_shared_tips(self):
        service = Tree()
        shared = service.shared_tips(["a", "b"], ["b", "c"])
        assert shared == ["b"]

    def test_shared_tips_no_overlap_exits(self, capsys):
        service = Tree()
        with pytest.raises(SystemExit) as exc:
            service.shared_tips(["a"], ["b"])
        assert exc.value.code == 2
        out, _ = capsys.readouterr()
        assert "no common tips" in out

    def test_prune_tree_using_taxa_list(self):
        service = Tree()
        tree = _TreeObj()
        res = service.prune_tree_using_taxa_list(tree, ["a", "b"])
        assert res is tree
        assert tree.pruned == ["a", "b"]

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
