import pytest
import sys
from math import isclose
from mock import patch, call
from pathlib import Path
from textwrap import dedent

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestTreenessOverRCV(object):
    @patch("builtins.print")
    def test_treeness_over_rcv0(self, mocked_print):
        expected_result = "0.35\t0.126\t0.36"
        testargs = [
            "phykit",
            "treeness_over_rcv",
            "-t",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-a",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_treeness_over_rcv1(self, mocked_print):
        expected_result = "2.1374\t0.7695\t0.36"
        testargs = [
            "phykit",
            "treeness_over_rcv",
            "-t",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            "-a",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_treeness_over_rcv_alias0(self, mocked_print):
        expected_result = "2.1374\t0.7695\t0.36"
        testargs = [
            "phykit",
            "toverr",
            "-t",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            "-a",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_treeness_over_rcv_alias1(self, mocked_print):
        expected_result = "2.1374\t0.7695\t0.36"
        testargs = [
            "phykit",
            "tor",
            "-t",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            "-a",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_treeness_over_rcv_incorrect_tree_path(self, mocked_print):
        testargs = [
            "phykit",
            "treeness_over_rcv",
            "-t",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tr",
            "-a",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_treeness_over_rcv_incorrect_alignment_path(self, mocked_print):
        testargs = [
            "phykit",
            "treeness_over_rcv",
            "-t",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            "-a",
            f"{here.parent.parent.parent}/sample_files/simple.f",
        ]
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2