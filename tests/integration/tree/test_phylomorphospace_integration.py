import json
import sys
from pathlib import Path

import pytest
from mock import patch

from phykit.phykit import Phykit


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
MULTI_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_multi_traits.tsv")


@pytest.mark.integration
class TestPhylomorphospaceIntegration:
    @patch("builtins.print")
    def test_basic_invocation(self, mocked_print, tmp_path):
        pytest.importorskip("matplotlib")
        plot_path = str(tmp_path / "phmo.png")
        testargs = [
            "phykit",
            "phylomorphospace",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "--trait-x", "body_mass",
            "--trait-y", "brain_size",
            "--plot-output", plot_path,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert Path(plot_path).exists()
        assert Path(plot_path).stat().st_size > 0

    @patch("builtins.print")
    def test_alias_phylomorpho(self, mocked_print, tmp_path):
        pytest.importorskip("matplotlib")
        plot_path = str(tmp_path / "phmo_alias.png")
        testargs = [
            "phykit",
            "phylomorpho",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "--trait-x", "body_mass",
            "--trait-y", "brain_size",
            "--plot-output", plot_path,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert Path(plot_path).exists()

    @patch("builtins.print")
    def test_alias_phmo(self, mocked_print, tmp_path):
        pytest.importorskip("matplotlib")
        plot_path = str(tmp_path / "phmo_short.png")
        testargs = [
            "phykit",
            "phmo",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "--trait-x", "body_mass",
            "--trait-y", "brain_size",
            "--plot-output", plot_path,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert Path(plot_path).exists()

    @patch("builtins.print")
    def test_with_color_by(self, mocked_print, tmp_path):
        pytest.importorskip("matplotlib")
        plot_path = str(tmp_path / "phmo_color.png")
        testargs = [
            "phykit",
            "phylomorphospace",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "--trait-x", "body_mass",
            "--trait-y", "brain_size",
            "--color-by", "longevity",
            "--plot-output", plot_path,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert Path(plot_path).exists()
        assert Path(plot_path).stat().st_size > 0

    @patch("builtins.print")
    def test_json_output(self, mocked_print, tmp_path):
        pytest.importorskip("matplotlib")
        plot_path = str(tmp_path / "phmo_json.png")
        testargs = [
            "phykit",
            "phylomorphospace",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "--trait-x", "body_mass",
            "--trait-y", "brain_size",
            "--plot-output", plot_path,
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["trait_x"] == "body_mass"
        assert payload["trait_y"] == "brain_size"
        assert "tip_data" in payload

    @patch("builtins.print")
    def test_auto_select_two_traits(self, mocked_print, tmp_path):
        pytest.importorskip("matplotlib")
        # Create a 2-trait file
        trait_file = tmp_path / "two_traits.tsv"
        trait_file.write_text(
            "taxon\ttrait_a\ttrait_b\n"
            "raccoon\t1.0\t2.0\n"
            "bear\t2.0\t3.0\n"
            "sea_lion\t3.0\t4.0\n"
            "seal\t4.0\t5.0\n"
            "monkey\t5.0\t6.0\n"
            "cat\t6.0\t7.0\n"
            "weasel\t7.0\t8.0\n"
            "dog\t8.0\t9.0\n"
        )
        plot_path = str(tmp_path / "phmo_auto.png")
        testargs = [
            "phykit",
            "phylomorphospace",
            "-t", TREE_SIMPLE,
            "-d", str(trait_file),
            "--plot-output", plot_path,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert Path(plot_path).exists()
