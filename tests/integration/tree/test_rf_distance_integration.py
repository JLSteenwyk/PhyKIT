import pytest
import sys
from math import isclose
from mock import patch, call
from pathlib import Path
from textwrap import dedent

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestTree(object):
    @patch("builtins.print")
    def test_rf_distance(self, mocked_print):
        expected_result = "8\t0.8"
        testargs = [
            "phykit",
            "rf_distance",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_other_topology.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_rf_distance_trees_with_different_tips(self, mocked_print):
        expected_result = "5\t0.8333"
        testargs = [
            "phykit",
            "rf_distance",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_other_topology_incomplete_taxon_representation.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_rf_distance_alias0(self, mocked_print):
        expected_result = "8\t0.8"
        testargs = [
            "phykit",
            "robinson_foulds_distance",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_other_topology.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_rf_distance_alias1(self, mocked_print):
        expected_result = "8\t0.8"
        testargs = [
            "phykit",
            "rf_dist",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_other_topology.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_rf_distance_alias2(self, mocked_print):
        expected_result = "8\t0.8"
        testargs = [
            "phykit",
            "rf",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_other_topology.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_rf_distance_incorrect_tree0_path(self, mocked_print):
        testargs = [
            "phykit",
            "rf",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tr",
            f"{here.parent.parent.parent}/sample_files/tree_simple_other_topology.tre",
        ]
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_rf_distance_incorrect_tree1_path(self, mocked_print):
        testargs = [
            "phykit",
            "rf",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_other_topology.tr",
        ]
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2