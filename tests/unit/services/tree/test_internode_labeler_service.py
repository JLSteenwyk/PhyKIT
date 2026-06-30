from argparse import Namespace
from io import StringIO
import subprocess
import sys

from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin
import pytest

from phykit.services.tree.internode_labeler import InternodeLabeler


def test_module_import_does_not_import_json_or_typing():
    code = """
import sys
import phykit.services.tree.internode_labeler as module
assert callable(module.print_json)
assert "json" not in sys.modules
assert "typing" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


@pytest.fixture
def args():
    return Namespace(tree="/some/path/to/file.tre", output=None)


class _Node:
    def __init__(self):
        self.confidence = None


class _Tree:
    def __init__(self):
        self.nonterms = [_Node(), _Node(), _Node()]

    def get_nonterminals(self):
        return self.nonterms


class TestInternodeLabeler:
    def test_init_sets_expected_attrs(self, args):
        service = InternodeLabeler(args)
        assert service.tree_file_path == args.tree
        assert service.output_file_path == "/some/path/to/file.tre.internode_labels.tre"
        assert service.json_output is False

    def test_process_args_honors_custom_output(self):
        args = Namespace(tree="/some/path/to/file.tre", output="/tmp/out.tre", json=True)
        service = InternodeLabeler(args)
        assert service.output_file_path == "/tmp/out.tre"
        assert service.json_output is True

    def test_add_labels_to_tree(self, args):
        service = InternodeLabeler(args)
        tree = _Tree()
        count = service.add_labels_to_tree(tree)
        assert count == 3
        assert [node.confidence for node in tree.get_nonterminals()] == [1, 2, 3]

    def test_add_labels_to_standard_tree_uses_direct_traversal(self, args, monkeypatch):
        service = InternodeLabeler(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")

        def fail_get_nonterminals(*_args, **_kwargs):
            raise AssertionError("generic nonterminal traversal should not be used")

        monkeypatch.setattr(TreeMixin, "get_nonterminals", fail_get_nonterminals)

        count = service.add_labels_to_tree(tree)

        assert count == 3
        assert tree.root.confidence == 1
        assert tree.root.clades[0].confidence == 2
        assert tree.root.clades[1].confidence == 3

    def test_add_labels_to_standard_tree_handles_mixed_child_counts(
        self,
        args,
        monkeypatch,
    ):
        service = InternodeLabeler(args)
        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);"),
            "newick",
        )

        def fail_get_nonterminals(*_args, **_kwargs):
            raise AssertionError("generic nonterminal traversal should not be used")

        monkeypatch.setattr(TreeMixin, "get_nonterminals", fail_get_nonterminals)

        count = service.add_labels_to_tree(tree)

        assert count == 3
        assert tree.root.confidence == 1
        assert tree.root.clades[1].confidence == 2
        assert tree.root.clades[2].confidence == 3

    def test_run_json_payload(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", output="/tmp/out.tre", json=True)
        service = InternodeLabeler(args)
        tree = _Tree()
        mocker.patch.object(InternodeLabeler, "read_tree_file", return_value=tree)
        mocked_write = mocker.patch.object(InternodeLabeler, "write_tree_file")
        mocked_json = mocker.patch("phykit.services.tree.internode_labeler.print_json")

        service.run()

        mocked_write.assert_called_once_with(tree, "/tmp/out.tre")
        payload = mocked_json.call_args.args[0]
        assert payload == {
            "input_tree": "/some/path/to/file.tre",
            "labeled_internodes": 3,
            "output_file": "/tmp/out.tre",
        }
