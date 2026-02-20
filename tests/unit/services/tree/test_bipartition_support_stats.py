import pytest
import sys
from argparse import Namespace
from mock import patch, call
import numpy as np

from phykit.phykit import Phykit
from phykit.services.tree.bipartition_support_stats import BipartitionSupportStats

@pytest.fixture
def args():
    kwargs = dict(tree="/some/path/to/file.tre", verbose=None, thresholds=None, json=False)
    return Namespace(**kwargs)


class TestBipartitionSupportStats(object):
    def test_init_sets_tree_file_path(self, args):
        t = BipartitionSupportStats(args)
        assert t.tree_file_path == args.tree
        assert t.output_file_path is None

    def test_read_file_reads_tree_file_path(self, mocker, args):
        # Mock the cached tree read method instead of Phylo.read
        mock_cached_read = mocker.patch("phykit.services.tree.base.Tree._cached_tree_read")
        mock_get_hash = mocker.patch("phykit.services.tree.base.Tree._get_file_hash", return_value="test_hash")

        t = BipartitionSupportStats(args)
        t.read_tree_file()

        # Verify the cached read was called with the correct parameters
        mock_get_hash.assert_called_with(args.tree)
        mock_cached_read.assert_called_with(args.tree, "newick", "test_hash")

    @patch("builtins.print")
    def test_bad_file_path(self, mocked_print):
        testargs = [
            "phykit",
            "bss",
            "some/file/path",
        ]

        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()

        assert pytest_wrapped_e.type is SystemExit
        mocked_print.assert_has_calls([
            call("Please check filename and pathing"),
        ])

    def test_parse_thresholds(self, args):
        args.thresholds = "70, 85.5,90"
        t = BipartitionSupportStats(args)
        assert t.thresholds == [70.0, 85.5, 90.0]

    def test_parse_thresholds_none(self, args):
        t = BipartitionSupportStats(args)
        assert t.parse_thresholds(None) == []

    def test_parse_thresholds_invalid_exits(self, args, capsys):
        t = BipartitionSupportStats(args)
        with pytest.raises(SystemExit) as exc:
            t.parse_thresholds("70,abc")
        assert exc.value.code == 2
        out, _ = capsys.readouterr()
        assert "must be numeric" in out

    def test_parse_thresholds_empty_exits(self, args, capsys):
        t = BipartitionSupportStats(args)
        with pytest.raises(SystemExit) as exc:
            t.parse_thresholds(" , ")
        assert exc.value.code == 2
        out, _ = capsys.readouterr()
        assert "Provide at least one numeric cutoff" in out

    def test_threshold_stats(self, args):
        args.thresholds = "90,100"
        t = BipartitionSupportStats(args)
        stats = t.calculate_threshold_stats([85.0, 85.0, 100.0, 100.0], t.thresholds)
        assert stats == [
            {"threshold": 90.0, "count_below": 2, "fraction_below": 0.5},
            {"threshold": 100.0, "count_below": 2, "fraction_below": 0.5},
        ]

    def test_threshold_stats_empty_bs_vals(self, args):
        args.thresholds = "70,90"
        t = BipartitionSupportStats(args)
        stats = t.calculate_threshold_stats([], t.thresholds)
        assert stats == [
            {"threshold": 70.0, "count_below": 0, "fraction_below": 0.0},
            {"threshold": 90.0, "count_below": 0, "fraction_below": 0.0},
        ]

    def test_to_builtin_converts_numpy_scalars(self, args):
        t = BipartitionSupportStats(args)
        value = {"a": np.int64(1), "b": [np.float64(1.5)]}
        converted = t._to_builtin(value)
        assert converted == {"a": 1, "b": [1.5]}

    def test_get_bipartition_support_vals(self, small_aspergillus_tree, args):
        t = BipartitionSupportStats(args)
        vals, names = t.get_bipartition_support_vals(small_aspergillus_tree)
        assert len(vals) == len(names)
        assert len(vals) > 0
        assert isinstance(names[0], list)

    def test_build_json_output_verbose(self, args):
        args.verbose = True
        t = BipartitionSupportStats(args)
        payload = t.build_json_output(
            bs_vals=[90.0, 100.0],
            term_names=[["A", "B"], ["C", "D"]],
            threshold_stats=[{"threshold": 95.0, "count_below": 1, "fraction_below": 0.5}],
        )
        assert payload["verbose"] is True
        assert payload["bipartitions"][0]["support"] == 90.0
        assert payload["thresholds"][0]["threshold"] == 95.0

    def test_build_json_output_non_verbose(self, args):
        args.verbose = False
        t = BipartitionSupportStats(args)
        payload = t.build_json_output(
            bs_vals=[90.0, 100.0],
            term_names=[["A", "B"], ["C", "D"]],
            threshold_stats=[],
        )
        assert payload["verbose"] is False
        assert "summary" in payload

    @patch("builtins.print")
    def test_json_output_mode(self, mocked_print, mocker, args):
        args.json = True
        args.verbose = False
        t = BipartitionSupportStats(args)
        mock_tree = mocker.Mock()
        mocker.patch.object(t, "read_tree_file", return_value=mock_tree)
        mocker.patch.object(
            t,
            "get_bipartition_support_vals",
            return_value=([85.0, 100.0], [["a", "b"], ["c", "d"]]),
        )
        dumped_json = '{"summary":{"mean":92.5},"verbose":false,"thresholds":[]}'
        mock_dumps = mocker.patch(
            "phykit.services.tree.bipartition_support_stats.json.dumps",
            return_value=dumped_json,
        )

        t.run()

        assert mock_dumps.called
        payload = mock_dumps.call_args.args[0]
        assert payload["verbose"] is False
        assert payload["summary"]["mean"] == 92.5
        mocked_print.assert_called_with(dumped_json)

    def test_json_output_mode_handles_broken_pipe(self, mocker, args):
        args.json = True
        args.verbose = True
        t = BipartitionSupportStats(args)
        mocker.patch.object(t, "read_tree_file", return_value=object())
        mocker.patch.object(t, "get_bipartition_support_vals", return_value=([85.0], [["a", "b"]]))
        mocker.patch.object(t, "calculate_threshold_stats", return_value=[])
        mocker.patch("phykit.services.tree.bipartition_support_stats.json.dumps", return_value='{"ok": true}')
        mocker.patch("builtins.print", side_effect=BrokenPipeError)
        t.run()

    def test_run_verbose_prints_bipartitions(self, mocker, args):
        args.verbose = True
        args.json = False
        t = BipartitionSupportStats(args)
        mocker.patch.object(t, "read_tree_file", return_value=object())
        mocker.patch.object(t, "get_bipartition_support_vals", return_value=([85.0], [["a", "b"]]))
        mocker.patch.object(t, "calculate_threshold_stats", return_value=[])
        mocked_print = mocker.patch("builtins.print")
        t.run()
        mocked_print.assert_called_with(85.0, "a;b")

    def test_run_nonverbose_prints_summary_and_thresholds(self, mocker, args):
        args.verbose = False
        args.json = False
        t = BipartitionSupportStats(args)
        mocker.patch.object(t, "read_tree_file", return_value=object())
        mocker.patch.object(t, "get_bipartition_support_vals", return_value=([85.0, 100.0], [["a"], ["b"]]))
        mocker.patch.object(
            t,
            "calculate_threshold_stats",
            return_value=[{"threshold": 90.0, "count_below": 1, "fraction_below": 0.5}],
        )
        mocked_summary = mocker.patch("phykit.services.tree.bipartition_support_stats.print_summary_statistics")
        mocked_print = mocker.patch("builtins.print")
        t.run()
        mocked_summary.assert_called_once()
        mocked_print.assert_called_with("below 90.0: 1 (50.0%)")


    # def test_calculate_bipartition_support_stats(self, small_aspergillus_tree, args):
    #     t = BipartitionSupportStats(args)
    #     bs_vals = t.get_bipartition_support_vals(small_aspergillus_tree)
    #     assert isinstance(stats["mean"], float)
    #     assert isinstance(stats["median"], (int or float))
    #     assert isinstance(stats["twenty_fifth"], float)
    #     assert isinstance(stats["seventy_fifth"], float)
    #     assert isinstance(stats["standard_deviation"], float)
    #     assert isinstance(stats["variance"], float)
    #     assert isclose(stats["mean"], 95.71428571428571, rel_tol=0.001)
    #     assert isclose(stats["median"], 100, rel_tol=0.001)
    #     assert isclose(stats["twenty_fifth"], 92.5, rel_tol=0.001)
    #     assert isclose(stats["seventy_fifth"], 100.0, rel_tol=0.001)
    #     assert isclose(stats["standard_deviation"], 7.319250547113999, rel_tol=0.001)
    #     assert isclose(stats["variance"], 53.57142857142857, rel_tol=0.001)
    #     assert isclose(stats["minimum"], 85, rel_tol=0.001)
    #     assert isclose(stats["maximum"], 100, rel_tol=0.001)
