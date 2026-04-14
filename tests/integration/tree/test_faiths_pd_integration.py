import pytest
import sys
import json
from mock import patch, call
from pathlib import Path

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestFaithsPD(object):
    @patch("builtins.print")
    def test_faiths_pd_include_root_default(self, mocked_print):
        expected = 62.2495
        testargs = [
            "phykit",
            "faiths_pd",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-t",
            f"{here.parent.parent.parent}/sample_files/faiths_pd_community_four.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected)]

    @patch("builtins.print")
    def test_faiths_pd_exclude_root_differs_when_mrca_below_root(self, mocked_print):
        expected = 26.0
        testargs = [
            "phykit",
            "faiths_pd",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-t",
            f"{here.parent.parent.parent}/sample_files/faiths_pd_community_mrca_deep.txt",
            "--exclude-root",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected)]

    @patch("builtins.print")
    def test_faiths_pd_include_root_when_mrca_below_root(self, mocked_print):
        expected = 26.846
        testargs = [
            "phykit",
            "faiths_pd",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-t",
            f"{here.parent.parent.parent}/sample_files/faiths_pd_community_mrca_deep.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected)]

    @patch("builtins.print")
    def test_faiths_pd_alias_fpd(self, mocked_print):
        expected = 62.2495
        testargs = [
            "phykit",
            "fpd",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-t",
            f"{here.parent.parent.parent}/sample_files/faiths_pd_community_four.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected)]

    @patch("builtins.print")
    def test_faiths_pd_json(self, mocked_print):
        testargs = [
            "phykit",
            "faiths_pd",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-t",
            f"{here.parent.parent.parent}/sample_files/faiths_pd_community_four.txt",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload == {"faiths_pd": 62.2495, "n_taxa": 4, "include_root": True}

    @patch("builtins.print")
    def test_faiths_pd_incorrect_tree_path(self, mocked_print):
        testargs = [
            "phykit",
            "faiths_pd",
            "/does/not/exist.tre",
            "-t",
            f"{here.parent.parent.parent}/sample_files/faiths_pd_community_four.txt",
        ]
        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as exc:
                Phykit()
        assert exc.value.code == 2

    @patch("builtins.print")
    def test_faiths_pd_incorrect_taxa_path(self, mocked_print):
        testargs = [
            "phykit",
            "faiths_pd",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-t",
            "/does/not/exist.txt",
        ]
        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as exc:
                Phykit()
        assert exc.value.code == 2
