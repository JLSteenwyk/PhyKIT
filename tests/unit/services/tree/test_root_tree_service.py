from argparse import Namespace
import subprocess
import sys

import pytest

from phykit.services.tree.root_tree import RootTree


@pytest.fixture
def args():
    return Namespace(
        tree="/some/path/to/file.tre",
        root="/some/path/to/outgroup.txt",
        output=None,
    )


class _Tree:
    pass


def test_tree_manipulation_modules_defer_heavy_imports():
    modules = [
        "phykit.services.tree.prune_tree",
        "phykit.services.tree.rename_tree_tips",
        "phykit.services.tree.root_tree",
        "phykit.services.tree.nearest_neighbor_interchange",
    ]
    code = f"""
import sys
modules = {modules!r}
for module_name in modules:
    module = __import__(module_name, fromlist=["*"])
    assert callable(module.print_json)

import phykit.services.tree.root_tree as root_module
assert hasattr(root_module.Phylo.BaseTree.Tree, "root_with_outgroup")
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "numpy" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_root_tree_import_defers_file_helper():
    code = """
import sys
import phykit.services.tree.root_tree as module
assert callable(module.read_single_column_file_to_list)
assert "phykit.helpers.files" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


class TestRootTree:
    def test_init_sets_expected_attrs(self, args):
        service = RootTree(args)
        assert service.tree_file_path == args.tree
        assert service.outgroup_taxa_file_path == args.root
        assert service.output_file_path == f"{args.tree}.rooted"
        assert service.json_output is False

    def test_process_args_honors_custom_output(self):
        args = Namespace(
            tree="/some/path/to/file.tre",
            root="/some/path/to/outgroup.txt",
            output="/tmp/rooted.tre",
            json=True,
        )
        service = RootTree(args)
        assert service.output_file_path == "/tmp/rooted.tre"
        assert service.json_output is True

    def test_run_roots_tree_and_writes_output(self, mocker):
        args = Namespace(
            tree="/some/path/to/file.tre",
            root="/some/path/to/outgroup.txt",
            output="/tmp/rooted.tre",
            json=False,
        )
        tree = _Tree()
        mocker.patch.object(RootTree, "read_tree_file", return_value=tree)
        mocker.patch(
            "phykit.services.tree.root_tree.read_single_column_file_to_list",
            return_value=["sea_lion", "seal"],
        )
        mocked_root = mocker.patch("phykit.services.tree.root_tree.Phylo.BaseTree.Tree.root_with_outgroup")
        mocked_write = mocker.patch.object(RootTree, "write_tree_file")

        service = RootTree(args)
        service.run()

        mocked_root.assert_called_once_with(tree, ["sea_lion", "seal"])
        mocked_write.assert_called_once_with(tree, "/tmp/rooted.tre")

    def test_run_json_outputs_payload(self, mocker):
        args = Namespace(
            tree="/some/path/to/file.tre",
            root="/some/path/to/outgroup.txt",
            output="/tmp/rooted.tre",
            json=True,
        )
        tree = _Tree()
        mocker.patch.object(RootTree, "read_tree_file", return_value=tree)
        mocker.patch(
            "phykit.services.tree.root_tree.read_single_column_file_to_list",
            return_value=["sea_lion", "seal"],
        )
        mocker.patch("phykit.services.tree.root_tree.Phylo.BaseTree.Tree.root_with_outgroup")
        mocker.patch.object(RootTree, "write_tree_file")
        mocked_json = mocker.patch("phykit.services.tree.root_tree.print_json")

        service = RootTree(args)
        service.run()

        payload = mocked_json.call_args.args[0]
        assert payload["input_tree"] == "/some/path/to/file.tre"
        assert payload["outgroup_taxa_file"] == "/some/path/to/outgroup.txt"
        assert payload["outgroup_taxa"] == ["sea_lion", "seal"]
        assert payload["output_file"] == "/tmp/rooted.tre"
