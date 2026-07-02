from argparse import Namespace
from io import StringIO
import subprocess
import sys

import pytest
from Bio import Phylo

from phykit.services.tree.tip_labels import TipLabels


@pytest.fixture
def args():
    return Namespace(tree="/some/path/to/file.tre")


class _Tip:
    def __init__(self, name):
        self.name = name


class _Tree:
    def __init__(self, tips):
        self._tips = tips

    def get_terminals(self):
        return [_Tip(tip) for tip in self._tips]


def test_module_import_does_not_import_json_or_heavy_tree_modules():
    code = """
import sys
import phykit.services.tree.tip_labels as module

assert callable(module.print_json)
assert "json" not in sys.modules
assert "typing" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "numpy" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


class TestTipLabels:
    def test_init_sets_expected_attrs(self, args):
        service = TipLabels(args)
        assert service.tree_file_path == args.tree
        assert service.json_output is False

    def test_run_prints_tip_names(self, mocker, args, capsys):
        mocker.patch.object(
            TipLabels,
            "read_tree_file_unmodified",
            return_value=_Tree(["a", "b", "c"]),
        )

        service = TipLabels(args)
        service.run()

        assert capsys.readouterr().out == "a\nb\nc\n"

    def test_run_uses_unmodified_tree_read(self, mocker, args):
        tree = _Tree(["a", "b"])
        service = TipLabels(args)
        read_tree = mocker.patch.object(
            service,
            "read_tree_file_unmodified",
            return_value=tree,
        )
        mocked_write = mocker.patch("phykit.services.tree.tip_labels.sys.stdout.write")

        service.run()

        read_tree.assert_called_once_with()
        mocked_write.assert_called_once_with("a\nb\n")

    def test_run_json_prints_structured_payload(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", json=True)
        mocker.patch.object(
            TipLabels,
            "read_tree_file_unmodified",
            return_value=_Tree(["a", "b"]),
        )
        mocked_json = mocker.patch("phykit.services.tree.tip_labels.print_json")

        service = TipLabels(args)
        service.run()

        payload = mocked_json.call_args.args[0]
        assert payload["rows"] == [{"taxon": "a"}, {"taxon": "b"}]
        assert payload["tips"] == ["a", "b"]

    def test_run_uses_fast_tip_names_for_parsed_tree(self, mocker):
        tree = Phylo.read(StringIO("((a:1,b:1):1,c:2);"), "newick")

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("get_terminals fallback should not be called")

        mocker.patch.object(TipLabels, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(tree, "get_terminals", side_effect=fail_get_terminals)
        mocked_write = mocker.patch("phykit.services.tree.tip_labels.sys.stdout.write")

        service = TipLabels(Namespace(tree="x.tre"))
        service.run()

        mocked_write.assert_called_once_with("a\nb\nc\n")

    def test_run_ignores_broken_pipe(self, mocker, args):
        mocker.patch.object(
            TipLabels,
            "read_tree_file_unmodified",
            return_value=_Tree(["a"]),
        )
        mocker.patch(
            "phykit.services.tree.tip_labels.sys.stdout.write",
            side_effect=BrokenPipeError,
        )

        service = TipLabels(args)
        service.run()
