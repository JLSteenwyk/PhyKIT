import pytest
import sys
from math import isclose
from mock import patch, call
from pathlib import Path
from textwrap import dedent

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestPTT(object):
    @patch("builtins.print")
    def test_polytomy_test0(self, mocked_print):
        testargs = [
            "phykit",
            "polytomy_test",
            "-t",
            f"{here.parent.parent.parent}/sample_files/test_trees.txt",
            "-g",
            f"{here.parent.parent.parent}/sample_files/test_trees_groups.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call("Gene Support Frequency Results"),
            call("=============================="),
            call("chi-squared: 20.0"),
            call("p-value: 4.5e-05"),
            call("total genes: 10"),
            call("0-1: 10"),
            call("0-2: 0"),
            call("1-2: 0"),
            # call("\nTriplet Results"),
            # call("==============="),
            # call("chi-squared: 320.0"),
            # call("p-value: 0.0"),
            # call("total triplets: 160"),
            # call("0-1: 160"),
            # call("0-2: 0"),
            # call("1-2: 0")
        ]

    @patch("builtins.print")
    def test_polytomy_test1(self, mocked_print):
        testargs = [
            "phykit",
            "polytomy_test",
            "-t",
            f"{here.parent.parent.parent}/sample_files/polyt_test_trees.txt",
            "-g",
            f"{here.parent.parent.parent}/sample_files/polyt_test_groups.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call("Gene Support Frequency Results"),
            call("=============================="),
            call("chi-squared: 2.0"),
            call("p-value: 0.367879"),
            call("total genes: 4"),
            call("0-1: 2"),
            call("0-2: 2"),
            call("1-2: 0"),
            # call("\nTriplet Results"),
            # call("==============="),
            # call("chi-squared: 2.0"),
            # call("p-value: 0.367879"),
            # call("total triplets: 4"),
            # call("0-1: 2"),
            # call("0-2: 2"),
            # call("1-2: 0")
        ]

    @patch("builtins.print")
    def test_polytomy_test_incorrect_group_path(self, mocked_print):
        testargs = [
            "phykit",
            "polytomy_test",
            "-t",
            f"{here.parent.parent.parent}/sample_files/polyt_test_trees.txt",
            "-g",
            f"{here.parent.parent.parent}/sample_files/polyt_test_groups.tx",
        ]
        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()

        assert pytest_wrapped_e.type == SystemExit
        mocked_print.assert_has_calls([
            call("Please check filename and pathing again."),
        ])
        

    @patch("builtins.print")
    def test_polytomy_test_incorrect_group_formatting(self, mocked_print):
        testargs = [
            "phykit",
            "polytomy_test",
            "-t",
            f"{here.parent.parent.parent}/sample_files/polyt_test_trees.txt",
            "-g",
            f"{here.parent.parent.parent}/sample_files/polyt_test_groups_three_tabs.txt",
        ]

        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()

        assert pytest_wrapped_e.type == SystemExit
        mocked_print.assert_has_calls([
            call("Please format the groups file (-g) as a four column tab-delimited file with column 1 being the name of the test"),
            call("col2: the tip names of one group (; separated)"),
            call("col3: the tip names of a second group (; separated)"),
            call("col4: the tip names of a third group (; separated)"),
            call("col5: the tip names of the outgroup taxa (; separated)"),
        ])
            

    @patch("builtins.print")
    def test_polytomy_test_incorrect_trees_path(self, mocked_print):
        testargs = [
            "phykit",
            "polytomy_test",
            "-t",
            f"{here.parent.parent.parent}/sample_files/polyt_test_trees.tx",
            "-g",
            f"{here.parent.parent.parent}/sample_files/polyt_test_groups.txt",
        ]
        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()

        assert pytest_wrapped_e.type == SystemExit
        mocked_print.assert_has_calls([
            call("Please check file name and pathing"),
        ])

    @patch("builtins.print")
    def test_polytomy_test_incorrect_path_within_trees_file(self, mocked_print):
        testargs = [
            "phykit",
            "polytomy_test",
            "-t",
            f"{here.parent.parent.parent}/sample_files/polyt_test_trees_bad_pathing.txt",
            "-g",
            f"{here.parent.parent.parent}/sample_files/polyt_test_groups.txt",
        ]
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_polytomy_test_incorrect_groups_file_format(self, mocked_print):
        testargs = [
            "phykit",
            "polytomy_test",
            "-t",
            f"{here.parent.parent.parent}/sample_files/polyt_test_trees.txt",
            "-g",
            f"{here.parent.parent.parent}/sample_files/polyt_test_incorrect_groups.txt",
        ]
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_polytomy_test_alias0(self, mocked_print):
        testargs = [
            "phykit",
            "polyt_test",
            "-t",
            f"{here.parent.parent.parent}/sample_files/polyt_test_trees.txt",
            "-g",
            f"{here.parent.parent.parent}/sample_files/polyt_test_groups.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call("Gene Support Frequency Results"),
            call("=============================="),
            call("chi-squared: 2.0"),
            call("p-value: 0.367879"),
            call("total genes: 4"),
            call("0-1: 2"),
            call("0-2: 2"),
            call("1-2: 0"),
            # call("\nTriplet Results"),
            # call("==============="),
            # call("chi-squared: 2.0"),
            # call("p-value: 0.367879"),
            # call("total triplets: 4"),
            # call("0-1: 2"),
            # call("0-2: 2"),
            # call("1-2: 0")
        ]

    @patch("builtins.print")
    def test_polytomy_test_alias1(self, mocked_print):
        testargs = [
            "phykit",
            "polyt",
            "-t",
            f"{here.parent.parent.parent}/sample_files/polyt_test_trees.txt",
            "-g",
            f"{here.parent.parent.parent}/sample_files/polyt_test_groups.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call("Gene Support Frequency Results"),
            call("=============================="),
            call("chi-squared: 2.0"),
            call("p-value: 0.367879"),
            call("total genes: 4"),
            call("0-1: 2"),
            call("0-2: 2"),
            call("1-2: 0"),
            # call("\nTriplet Results"),
            # call("==============="),
            # call("chi-squared: 2.0"),
            # call("p-value: 0.367879"),
            # call("total triplets: 4"),
            # call("0-1: 2"),
            # call("0-2: 2"),
            # call("1-2: 0")
        ]

    @patch("builtins.print")
    def test_polytomy_test_alias2(self, mocked_print):
        testargs = [
            "phykit",
            "ptt",
            "-t",
            f"{here.parent.parent.parent}/sample_files/polyt_test_trees.txt",
            "-g",
            f"{here.parent.parent.parent}/sample_files/polyt_test_groups.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call("Gene Support Frequency Results"),
            call("=============================="),
            call("chi-squared: 2.0"),
            call("p-value: 0.367879"),
            call("total genes: 4"),
            call("0-1: 2"),
            call("0-2: 2"),
            call("1-2: 0"),
            # call("\nTriplet Results"),
            # call("==============="),
            # call("chi-squared: 2.0"),
            # call("p-value: 0.367879"),
            # call("total triplets: 4"),
            # call("0-1: 2"),
            # call("0-2: 2"),
            # call("1-2: 0")
        ]
