import pytest
from argparse import Namespace
from Bio import Phylo
from math import isclose

from phykit.services.tree.internode_labeler import InternodeLabeler


class TestInternodeLabeler(object):
    def test_init_sets_tree_file_path(self, args):
        internode_labeler = InternodeLabeler(args)
        assert internode_labeler.tree_file_path == args.tree

    def test_init_sets_output_file_path(self, args):
        internode_labeler = InternodeLabeler(args)
        assert internode_labeler.output_file_path == f"{args.tree}.internodeLabels.tree"

    def test_read_file_reads_tree_file_path(self, mocker, args):
        mock_read = mocker.patch("phykit.services.tree.base.Phylo.read")
        internode_labeler = InternodeLabeler(args)
        internode_labeler.read_tree_file()
        mock_read.assert_called_with(args.tree, "newick")

    def test_add_labels_to_tree(self, tree_simple, args):
        internode_labeler = InternodeLabeler(args)
        for node in tree_simple.get_nonterminals():
            assert node.confidence is None

        labeled_tree = internode_labeler.add_labels_to_tree(tree_simple)

        for node in labeled_tree.get_nonterminals():
            assert isinstance(node.confidence, int)
