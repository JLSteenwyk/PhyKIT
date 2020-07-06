import pytest
import sys
from math import isclose
from mock import patch, call
from pathlib import Path
from textwrap import dedent

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestDVMC(object):
    @patch("builtins.print")
    def test_dvmc(self, mocked_print):
        expected_result = 42.8016
        testargs = [
            "phykit",
            "dvmc",
            "-t",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-r",
            f"{here.parent.parent.parent}/sample_files/tree_simple.outgroup.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_dvmc_alias(self, mocked_print):
        expected_result = 42.8016
        testargs = [
            "phykit",
            "degree_of_violation_of_a_molecular_clock",
            "-t",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-r",
            f"{here.parent.parent.parent}/sample_files/tree_simple.outgroup.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_dvmc_incorrect_outgroup_file(self, mocked_print):
        testargs = [
            "phykit",
            "dvmc",
            "-t",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-r",
            f"{here.parent.parent.parent}/sample_files/tree_simple.outgroup.tx",
        ]
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_dvmc_incorrect_tree_file(self, mocked_print):
        testargs = [
            "phykit",
            "dvmc",
            "-t",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tr",
            "-r",
            f"{here.parent.parent.parent}/sample_files/tree_simple.outgroup.txt",
        ]
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2