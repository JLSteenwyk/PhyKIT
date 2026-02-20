from argparse import Namespace

import pytest

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


class TestTipLabels:
    def test_init_sets_expected_attrs(self, args):
        service = TipLabels(args)
        assert service.tree_file_path == args.tree
        assert service.json_output is False

    def test_run_prints_tip_names(self, mocker, args):
        mocker.patch.object(TipLabels, "read_tree_file", return_value=_Tree(["a", "b", "c"]))
        mocked_print = mocker.patch("builtins.print")

        service = TipLabels(args)
        service.run()

        mocked_print.assert_called_once_with("a\nb\nc")

    def test_run_json_prints_structured_payload(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", json=True)
        mocker.patch.object(TipLabels, "read_tree_file", return_value=_Tree(["a", "b"]))
        mocked_json = mocker.patch("phykit.services.tree.tip_labels.print_json")

        service = TipLabels(args)
        service.run()

        payload = mocked_json.call_args.args[0]
        assert payload["rows"] == [{"taxon": "a"}, {"taxon": "b"}]
        assert payload["tips"] == ["a", "b"]

    def test_run_ignores_broken_pipe(self, mocker, args):
        mocker.patch.object(TipLabels, "read_tree_file", return_value=_Tree(["a"]))
        mocker.patch("builtins.print", side_effect=BrokenPipeError)

        service = TipLabels(args)
        service.run()
