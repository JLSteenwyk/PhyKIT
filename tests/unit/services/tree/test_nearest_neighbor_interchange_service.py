from argparse import Namespace
from io import StringIO

import pytest
from Bio import Phylo

from phykit.services.tree.nearest_neighbor_interchange import (
    NearestNeighborInterchange,
)


@pytest.fixture
def args():
    return Namespace(tree="/some/path/to/file.tre", output=None)


class TestNearestNeighborInterchange:
    def test_init_sets_expected_attrs(self, args):
        service = NearestNeighborInterchange(args)
        assert service.tree_file_path == args.tree
        assert service.output_file_path == f"{args.tree}.nnis"
        assert service.json_output is False

    def test_process_args_honors_custom_output(self):
        args = Namespace(tree="/some/path/to/file.tre", output="/tmp/nnis.tre", json=True)
        service = NearestNeighborInterchange(args)
        assert service.output_file_path == "/tmp/nnis.tre"
        assert service.json_output is True

    def test_fast_tree_copy_returns_distinct_object(self, args):
        service = NearestNeighborInterchange(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        copied = service._fast_tree_copy(tree)
        assert copied is not tree
        assert copied.format("newick") == tree.format("newick")

    def test_get_neighbors_returns_list(self, args):
        service = NearestNeighborInterchange(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        neighbors = service.get_neighbors(tree)
        assert isinstance(neighbors, list)
        assert len(neighbors) >= 1

    def test_get_neighbors_exercises_right_child_branch(self, args):
        service = NearestNeighborInterchange(args)
        # This topology has a non-root internal clade as parent's right child,
        # which exercises the alternate branch in get_neighbors.
        tree = Phylo.read(StringIO("((A:1,(B:1,C:1):1):1,(D:1,E:1):1);"), "newick")
        neighbors = service.get_neighbors(tree)
        assert len(neighbors) >= 1
        assert all(hasattr(n, "root") for n in neighbors)

    def test_get_neighbors_exercises_left_child_branch(self, args):
        service = NearestNeighborInterchange(args)
        # This topology has a non-root internal clade as parent's left child.
        tree = Phylo.read(StringIO("((((A:1,B:1):1,C:1):1,D:1):1,(E:1,F:1):1);"), "newick")
        neighbors = service.get_neighbors(tree)
        assert len(neighbors) >= 1

    def test_run_json_payload(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", output="/tmp/nnis.tre", json=True)
        service = NearestNeighborInterchange(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        mocker.patch.object(NearestNeighborInterchange, "read_tree_file", return_value=tree)
        mocker.patch.object(NearestNeighborInterchange, "get_neighbors", return_value=[tree, tree])
        mocked_write = mocker.patch("phykit.services.tree.nearest_neighbor_interchange.Phylo.write")
        mocked_json = mocker.patch("phykit.services.tree.nearest_neighbor_interchange.print_json")

        service.run()

        mocked_write.assert_called_once()
        payload = mocked_json.call_args.args[0]
        assert payload == {
            "input_tree": "/some/path/to/file.tre",
            "total_trees": 3,
            "nni_neighbors": 2,
            "output_file": "/tmp/nnis.tre",
        }
