from argparse import Namespace

import pytest

from phykit.services.tree.monophyly_check import MonophylyCheck


@pytest.fixture
def args():
    return Namespace(tree="/some/path/to/file.tre", list_of_taxa="/some/path/to/taxa.txt")


class _Node:
    def __init__(self, confidence):
        self.confidence = confidence


class _Clade:
    def get_nonterminals(self):
        return [_Node(90), _Node(80), _Node(None)]


class TestMonophylyCheck:
    def test_init_sets_expected_attrs(self, args):
        service = MonophylyCheck(args)
        assert service.tree_file_path == args.tree
        assert service.list_of_taxa == args.list_of_taxa
        assert service.json_output is False

    def test_get_bootstrap_statistics(self, args):
        service = MonophylyCheck(args)
        stats = service.get_bootstrap_statistics(_Clade())
        assert round(stats["mean"], 4) == 85.0
        assert stats["maximum"] == 90
        assert stats["minimum"] == 80

    def test_populate_res_arr_monophyletic(self, args):
        service = MonophylyCheck(args)
        stats = {"mean": 90.0, "maximum": 100.0, "minimum": 80.0, "standard_deviation": 10.0}
        res = service.populate_res_arr([], stats, [])
        assert res[0][0] == "monophyletic"
        assert res[0][1:5] == [90.0, 100.0, 80.0, 10.0]

    def test_print_results_json(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", list_of_taxa="/some/path/to/taxa.txt", json=True)
        service = MonophylyCheck(args)
        mocked_json = mocker.patch("phykit.services.tree.monophyly_check.print_json")
        service.print_results(
            [["not_monophyletic", 95.0, 100.0, 85.0, 7.0, ["z", "a"]]]
        )
        payload = mocked_json.call_args.args[0]
        assert payload["rows"] == payload["results"]
        assert payload["rows"][0]["status"] == "not_monophyletic"
        assert payload["rows"][0]["offending_taxa"] == ["a", "z"]

    def test_run_insufficient_taxa_exits(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", list_of_taxa="/some/path/to/taxa.txt")
        service = MonophylyCheck(args)
        mocker.patch.object(MonophylyCheck, "read_tree_file", return_value=object())
        mocker.patch(
            "phykit.services.tree.monophyly_check.read_single_column_file_to_list",
            return_value=["only_one_taxon"],
        )
        mocker.patch.object(MonophylyCheck, "get_tip_names_from_tree", return_value=["only_one_taxon"])
        with pytest.raises(SystemExit) as excinfo:
            service.run()
        assert excinfo.value.code == 2
