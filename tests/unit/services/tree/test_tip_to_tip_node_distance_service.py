from argparse import Namespace

import pytest

from phykit.services.tree.tip_to_tip_node_distance import TipToTipNodeDistance


@pytest.fixture
def args():
    return Namespace(tree_zero="/some/path/to/file.tre", tip_1="a", tip_2="b")


class _Tree:
    pass


class TestTipToTipNodeDistance:
    def test_init_sets_expected_attrs(self, args):
        service = TipToTipNodeDistance(args)
        assert service.tree_file_path == args.tree_zero
        assert service.tip_1 == "a"
        assert service.tip_2 == "b"
        assert service.json_output is False

    def test_check_leaves_exits_when_first_missing(self, mocker, args):
        service = TipToTipNodeDistance(args)
        mocker.patch(
            "phykit.services.tree.tip_to_tip_node_distance.TreeMixin.find_any",
            side_effect=[None, object()],
        )
        with pytest.raises(SystemExit) as excinfo:
            service.check_leaves(_Tree(), "missing", "b")
        assert excinfo.value.code == 2

    def test_run_text_output(self, mocker):
        args = Namespace(tree_zero="/some/path/to/file.tre", tip_1="a", tip_2="b", json=False)
        service = TipToTipNodeDistance(args)
        mocker.patch.object(TipToTipNodeDistance, "read_tree_file", return_value=_Tree())
        mocker.patch.object(TipToTipNodeDistance, "check_leaves")
        mocker.patch("phykit.services.tree.tip_to_tip_node_distance.TreeMixin.trace", return_value=[1, 2, 3, 4])
        mocked_print = mocker.patch("builtins.print")
        service.run()
        mocked_print.assert_called_once_with(4)

    def test_run_json_output(self, mocker):
        args = Namespace(tree_zero="/some/path/to/file.tre", tip_1="a", tip_2="b", json=True)
        service = TipToTipNodeDistance(args)
        mocker.patch.object(TipToTipNodeDistance, "read_tree_file", return_value=_Tree())
        mocker.patch.object(TipToTipNodeDistance, "check_leaves")
        mocker.patch("phykit.services.tree.tip_to_tip_node_distance.TreeMixin.trace", return_value=[1, 2, 3])
        mocked_json = mocker.patch("phykit.services.tree.tip_to_tip_node_distance.print_json")
        service.run()
        payload = mocked_json.call_args.args[0]
        assert payload == {
            "taxon_a": "a",
            "taxon_b": "b",
            "tip_to_tip_node_distance": 3,
        }
